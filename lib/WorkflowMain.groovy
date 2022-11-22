// This file is based on the nf-core/tools pipeline-template.
// Changes to this file must be propagated via wf-template.

class WorkflowMain {

    // Citation string for pipeline
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n"
    }

    // Generate help string
    public static String help(workflow, params, log) {
        String line_sep = ' \\ \n\t'
        String command_example = params.wf.example_cmd.join(line_sep)
        String command = 'nextflow run ' + workflow.manifest.name + line_sep + command_example
        String help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        return help_string
    }

    // Generate parameter summary log string
    public static String paramsSummaryLog(workflow, params, log) {
        String workflow_version = NfcoreTemplate.version(workflow)
        String summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        summary_log += "\nThis is ${workflow.manifest.name} ${workflow_version}.\n"
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    // Validate parameters and print summary to screen
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Explode on conda
        // conda.enabled seems to be backward compatible but wrap this
        // in a generic catch just in case
        try {
            if (workflow.session.config.conda.enabled) {
                log.error "Sorry, this workflow is not compatible with Conda, please use -profile standard (Docker) or -profile singularity."
                System.exit(1)
            }
        } catch(Exception e) {}

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)
    }
}
