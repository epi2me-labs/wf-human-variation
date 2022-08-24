import static groovy.json.JsonOutput.toJson
import groovy.json.JsonBuilder
import groovy.json.JsonSlurper


class Pinguscript {
      public static String ping_post(workflow, message, error_message, out_dir, params) {
        def msgId = UUID.randomUUID().toString();
        def hosthash = null;
        try { 
            hosthash = InetAddress.getLocalHost().getHostName().md5()
        } catch(Exception e) {
            hosthash = "Unavailable"
        }
        def opsys = System.properties['os.name'].toLowerCase()
        if (System.properties['os.version'].toLowerCase().contains("wsl")){
            opsys = "WSL"
        }
        def workflow_name = "$workflow.manifest.name"
        def session = "$workflow.sessionId" 
        def errorMessage = "$error_message"
        def profile = "$workflow.profile"
        def filename =  "$out_dir/params.json"
        File fileb = new File(filename)
        def any_other_data = [:]
        if (fileb.exists() && "$message" != "start") {
            def jsonSlurper = new JsonSlurper()
            any_other_data = jsonSlurper.parse(fileb)
        } 
        def meta_json = new JsonBuilder()
        def agent = "$params.wf.agent"
        def meta = meta_json "error": errorMessage.toString(), "profile": profile.toString(),
            "agent": agent.toString()
        meta+=any_other_data
        def ping_version = '2.0.0'
        def tracking_json = new JsonBuilder()
        def tracking_id = tracking_json "msg_id": msgId, "version": ping_version
        def data_json = new JsonBuilder()
        def data = data_json "source": "workflow", "workflow": workflow_name.toString(),
                "message": message, "meta": meta
        def body_json = new JsonBuilder()
        def root = body_json "tracking_id": tracking_id,  "hostname": hosthash.toString(), "os": opsys.toString(),
                "session": session.toString(), "data": data
        String postResult
        ((HttpURLConnection)new URL('https://ping.oxfordnanoportal.com/epilaby').openConnection()).with({
            requestMethod = 'POST'
            doOutput = true
            setRequestProperty('Content-Type', 'application/json') 
            setRequestProperty('accept', 'application/json') 
            outputStream.withPrintWriter({printWriter ->
                 printWriter.write(body_json.toString())
            })
            postResult = inputStream.text
        })
        return (postResult)
    }
}