import static groovy.json.JsonOutput.toJson
import groovy.json.JsonBuilder
import groovy.json.JsonSlurper


class Pinguscript {

    // Send a ping for the start of a workflow
    public static void ping_start(nextflow, workflow, params) {
        wf_ping(nextflow, workflow, "start", null, params)
    }
    // Send a ping for a completed workflow (successful or otherwise)
    public static void ping_complete(nextflow, workflow, params) {
        wf_ping(nextflow, workflow, "end", null, params)
    }
    // Send a ping for a workflow error
    public static void ping_error(nextflow, workflow, params) {
        def error_message = workflow.errorMessage
        wf_ping(nextflow, workflow, "error", error_message, params)
    }
    // Shared handler to construct a ping JSON and send it
    private static String wf_ping(nextflow, workflow, event, error_message, params) {
        if (params.disable_ping) {
            return "{}"
        }
        def body_json = make_wf_ping(nextflow, workflow, event, error_message, params)
        send_ping_post("epilaby", body_json)
    }

    // Helper to removing keys from a map
    private static clean_meta(meta, keys_to_remove) {
        for (key in keys_to_remove) {
            if (meta.containsKey(key)) {
                meta.remove(key)
            }
        }
    }

    // Helper for fetching a key from the params map
    // seems pointless but you just know someone is going to end up writing meta.this ? meta.that
    private static get_meta(meta, key) {
        (meta.containsKey(key) && meta[key]) ? meta[key].toString() : null
    }

    // Construct workflow ping JSON
    private static String make_wf_ping(nextflow, workflow, event, error_message, params) {
        // cheeky deepcopy using json
        String paramsJSON = new JsonBuilder(params).toPrettyString()
        def params_data = new JsonSlurper().parseText(paramsJSON)

        // hostname
        def host = null
        try {
            host = InetAddress.getLocalHost().getHostName()
        }
        catch(Exception e) {}

        // OS
        // TODO check version on WSL
        def opsys = System.properties['os.name'].toLowerCase()
        def opver = System.properties['os.version']
        if (opver.toLowerCase().contains("wsl")){
            opsys = "wsl"
        }

        // placeholder for any future okta business
        // for now we'll use the guest_<ulid> sent to wf.epi2me_user
        def user = get_meta(params.wf, "epi2me_user")

        // drop cruft to save some precious bytes
        // affects the deep copy rather than original params
        clean_meta(params_data, [
            "schema_ignore_params",
        ])
        def ingress_ids = []
        if (params_data.containsKey("wf")) {
            ingress_ids = params_data.wf["ingress.run_ids"] ?: []
            clean_meta(params_data.wf, [
                "agent", // we send this later
                "epi2me_instance", // we send this later
                "epi2me_user", // we send this later
                "example_cmd",
                "ingress.run_ids", // we will send this elsewhere
            ])
        }

        // try and get runtime information
        def cpus = null
        try {
            cpus = Runtime.getRuntime().availableProcessors()
        }
        catch(Exception e) {}

        def workflow_success = null
        def workflow_exitcode = null
        if (event != "start") {
            workflow_success = workflow.success
            workflow_exitcode = workflow.exitStatus
        }

        /// build message
        def body_json = new JsonBuilder()
        body_json \
            "tracking_id": [
                "msg_id": UUID.randomUUID().toString(),
                "version": "3.0.0"
            ],
            "source": "workflow",
            "event": event,
            "params": params_data,
            // data will be null on start events, as ingress has not run
            "data": event != "start" ? [run_ids: ingress_ids] : null,
            "workflow": [
                "name": workflow.manifest.name,
                "version": workflow.manifest.version, // could use NfcoreTemplate.version(workflow)
                "run_name": workflow.runName, // required to disambiguate sessions
                "session": workflow.sessionId,
                "profile": workflow.profile,
                "resume": workflow.resume,
                "error": error_message, // null if no error
                "success": workflow_success,
                "exitcode": workflow_exitcode,
            ],
            "env": [
                "user": user, // placeholder for any future okta
                "hostname": host,
                "os": [
                    "name": opsys,
                    "version": opver
                ],
                "resource": [
                    "cpus": cpus,
                    "memory": null, // placeholder, no point asking via Runtime as it will just give us the Xmx size
                ],
                "agent": get_meta(params.wf, "agent"), // access via original params
                "epi2me": [
                    "instance": get_meta(params.wf, "epi2me_instance"),
                    "user": user,
                ],
                "nextflow": [
                    "version": nextflow.version.toString(),
                    "version_compat": nextflow.version.matches(workflow.manifest.nextflowVersion)
                ]
            ]
        return body_json
    }

    // Send a JSON payload to a given endpoint
    private static String send_ping_post(endpoint, body_json) {
        // Attempt to send payload and absorb any possible Exception gracefully
        String postResult
        boolean raise_exception = false
        try {
            ((HttpURLConnection)new URL("https://ping.oxfordnanoportal.com/${endpoint}").openConnection()).with({
                requestMethod = 'POST'
                doOutput = true
                setConnectTimeout(5000)
                setReadTimeout(10000)
                setRequestProperty('Content-Type', 'application/json')
                setRequestProperty('accept', 'application/json')
                outputStream.withPrintWriter({printWriter ->
                    printWriter.write(body_json.toString())
                })

                // Rethrow exceptions that imply we're not using this endpoint properly
                if(responseCode >= 400 && agent.toString() == "cw-ci") {
                    raise_exception = true
                }
                // Accessing inputStream.text will raise an Exception for failed requests
                postResult = inputStream.text
            })
        }
        catch(Exception e) {
            if(raise_exception) { throw e }
        }
        return (postResult)
    }
}
