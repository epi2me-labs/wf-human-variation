// This file is based on the nf-core/tools pipeline-template.
// Changes to this file must be propagated via wf-template.

class WorkflowMain {

    // Citation string for pipeline
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n"
    }

    // Print help to screen
    public static String help(workflow, params, log) {
        String line_sep = ' \\ \n\t'
        def command_example = params.wf.example_cmd.join(line_sep)
        def command = 'nextflow run ' + workflow.manifest.name + line_sep + command_example
        def help_string = ''
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        return help_string
    }

    // Print parameter summary log to screen
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        return summary_log
    }

    // Validate parameters and print summary to screen
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)
    }
}
