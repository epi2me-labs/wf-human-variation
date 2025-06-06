//
// This file holds several functions used within the nf-core pipeline template.
//

// MIT License
// 
// Copyright (c) 2018 nf-core
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


import org.yaml.snakeyaml.Yaml

class NfcoreTemplate {

    //
    // Check AWS Batch related parameters have been specified correctly
    //
    public static void awsBatch(workflow, params) {
        if (workflow.profile.contains('awsbatch')) {
            // Check params.awsqueue and params.awsregion have been set if running on AWSBatch
            assert (params.awsqueue && params.awsregion) : "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
            // Check outdir paths to be S3 buckets if running on AWSBatch
            assert params.outdir.startsWith('s3:')       : "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
        }
    }

    //
    // Check params.hostnames
    //
    public static void hostName(workflow, params, log) {
        Map colors = logColours(params.monochrome_logs)
        if (params.hostnames) {
            try {
                def hostname = "hostname".execute().text.trim()
                params.hostnames.each { prof, hnames ->
                    hnames.each { hname ->
                        if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                            log.info "=${colors.yellow}====================================================${colors.reset}=\n" +
                                "${colors.yellow}WARN: You are running with `-profile $workflow.profile`\n" +
                                "      but your machine hostname is ${colors.white}'$hostname'${colors.reset}.\n" +
                                "      ${colors.yellow_bold}Please use `-profile $prof${colors.reset}`\n" +
                                "=${colors.yellow}====================================================${colors.reset}="
                        }
                    }
                }
            } catch (Exception e) {
                log.warn "[$workflow.manifest.name] Could not determine 'hostname' - skipping check. Reason: ${e.message}."
            }
        }
    }

    //
    // Generate version string
    //
    public static String version(workflow) {
        String version_string = ""

        if (workflow.manifest.version) {
            def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
            version_string += "${prefix_v}${workflow.manifest.version}"
        }

        if (workflow.commitId) {
            def git_shortsha = workflow.commitId.substring(0, 7)
            version_string += "-g${git_shortsha}"
        }

        return version_string
    }

    //
    // Construct and send completion email
    //
    public static void email(workflow, params, summary_params, projectDir, log, multiqc_report=[], fail_mapped_reads=[:]) {

        // Set up the e-mail variables
        def subject = "[$workflow.manifest.name] Successful: $workflow.runName"
        if (fail_mapped_reads.size() > 0) {
            subject = "[$workflow.manifest.name] Partially successful (${fail_mapped_reads.size()} skipped): $workflow.runName"
        }
        if (!workflow.success) {
            subject = "[$workflow.manifest.name] FAILED: $workflow.runName"
        }

        def summary = [:]
        for (group in summary_params.keySet()) {
            summary << summary_params[group]
        }

        def misc_fields = [:]
        misc_fields['Date Started']              = workflow.start
        misc_fields['Date Completed']            = workflow.complete
        misc_fields['Pipeline script file path'] = workflow.scriptFile
        misc_fields['Pipeline script hash ID']   = workflow.scriptId
        if (workflow.repository) misc_fields['Pipeline repository Git URL']    = workflow.repository
        if (workflow.commitId)   misc_fields['Pipeline repository Git Commit'] = workflow.commitId
        if (workflow.revision)   misc_fields['Pipeline Git branch/tag']        = workflow.revision
        misc_fields['Nextflow Version']           = workflow.nextflow.version
        misc_fields['Nextflow Build']             = workflow.nextflow.build
        misc_fields['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

        def email_fields = [:]
        email_fields['version']           = NfcoreTemplate.version(workflow)
        email_fields['runName']           = workflow.runName
        email_fields['success']           = workflow.success
        email_fields['dateComplete']      = workflow.complete
        email_fields['duration']          = workflow.duration
        email_fields['exitStatus']        = workflow.exitStatus
        email_fields['errorMessage']      = (workflow.errorMessage ?: 'None')
        email_fields['errorReport']       = (workflow.errorReport ?: 'None')
        email_fields['commandLine']       = workflow.commandLine
        email_fields['projectDir']        = workflow.projectDir
        email_fields['summary']           = summary << misc_fields
        email_fields['fail_mapped_reads'] = fail_mapped_reads.keySet()
        email_fields['min_mapped_reads']  = params.min_mapped_reads

        // On success try attach the multiqc report
        def mqc_report = null
        try {
            if (workflow.success && !params.skip_multiqc) {
                mqc_report = multiqc_report.getVal()
                if (mqc_report.getClass() == ArrayList && mqc_report.size() >= 1) {
                    if (mqc_report.size() > 1) {
                        log.warn "[$workflow.manifest.name] Found multiple reports from process 'MULTIQC', will use only one"
                    }
                    mqc_report = mqc_report[0]
                }
            }
        } catch (all) {
            if (multiqc_report) {
                log.warn "[$workflow.manifest.name] Could not attach MultiQC report to summary email"
            }
        }

        // Check if we are only sending emails on failure
        def email_address = params.email
        if (!params.email && params.email_on_fail && !workflow.success) {
            email_address = params.email_on_fail
        }

        // Render the TXT template
        def engine       = new groovy.text.GStringTemplateEngine()
        def tf           = new File("$projectDir/assets/email_template.txt")
        def txt_template = engine.createTemplate(tf).make(email_fields)
        def email_txt    = txt_template.toString()

        // Render the HTML template
        def hf            = new File("$projectDir/assets/email_template.html")
        def html_template = engine.createTemplate(hf).make(email_fields)
        def email_html    = html_template.toString()

        // Render the sendmail template
        def max_multiqc_email_size = params.max_multiqc_email_size as nextflow.util.MemoryUnit
        def smail_fields           = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: max_multiqc_email_size.toBytes() ]
        def sf                     = new File("$projectDir/assets/sendmail_template.txt")
        def sendmail_template      = engine.createTemplate(sf).make(smail_fields)
        def sendmail_html          = sendmail_template.toString()

        // Send the HTML e-mail
        Map colors = logColours(params.monochrome_logs)
        if (email_address) {
            try {
                if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
                // Try to send HTML e-mail using sendmail
                [ 'sendmail', '-t' ].execute() << sendmail_html
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Sent summary e-mail to $email_address (sendmail)-"
            } catch (all) {
                // Catch failures and try with plaintext
                def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
                if ( mqc_report.size() <= max_multiqc_email_size.toBytes() ) {
                    mail_cmd += [ '-A', mqc_report ]
                }
                mail_cmd.execute() << email_html
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Sent summary e-mail to $email_address (mail)-"
            }
        }

        // Write summary e-mail HTML to a file
        def output_d = new File("${params.outdir}/pipeline_info/")
        if (!output_d.exists()) {
            output_d.mkdirs()
        }
        def output_hf = new File(output_d, "pipeline_report.html")
        output_hf.withWriter { w -> w << email_html }
        def output_tf = new File(output_d, "pipeline_report.txt")
        output_tf.withWriter { w -> w << email_txt }
    }

    //
    // Print pipeline summary on completion
    //
    public static void summary(workflow, params, log, fail_mapped_reads=[:], pass_mapped_reads=[:]) {
        Map colors = logColours(params.monochrome_logs)

        if (pass_mapped_reads.size() > 0) {
            def idx = 0
            def samp_aln = ''
            def total_aln_count = pass_mapped_reads.size() + fail_mapped_reads.size()
            for (samp in pass_mapped_reads) {
                samp_aln += "    ${samp.value}: ${samp.key}\n"
                idx += 1
                if (idx > 5) {
                    samp_aln += "    ..see pipeline reports for full list\n"
                    break;
                }
            }
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} ${pass_mapped_reads.size()}/$total_aln_count samples passed Bowtie2 ${params.min_mapped_reads} mapped read threshold:\n${samp_aln}${colors.reset}-"
        }
        if (fail_mapped_reads.size() > 0) {
            def samp_aln = ''
            for (samp in fail_mapped_reads) {
                samp_aln += "    ${samp.value}: ${samp.key}\n"
            }
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} ${fail_mapped_reads.size()} samples skipped since they failed Bowtie2 ${params.min_mapped_reads} mapped read threshold:\n${samp_aln}${colors.reset}-"
        }

        if (workflow.success) {
            if (workflow.stats.ignoredCount == 0) {
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}-"
            } else {
                log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed successfully, but with errored process(es) ${colors.reset}-"
            }
        } else {
            hostName(workflow, params, log)
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}-"
        }
    }

    //
    // ANSII Colours used for terminal logging
    //
    public static Map logColours(Boolean monochrome_logs) {
        Map colorcodes = [:]

        // Reset / Meta
        colorcodes['reset']      = monochrome_logs ? '' : "\033[0m"
        colorcodes['bold']       = monochrome_logs ? '' : "\033[1m"
        colorcodes['dim']        = monochrome_logs ? '' : "\033[2m"
        colorcodes['underlined'] = monochrome_logs ? '' : "\033[4m"
        colorcodes['blink']      = monochrome_logs ? '' : "\033[5m"
        colorcodes['reverse']    = monochrome_logs ? '' : "\033[7m"
        colorcodes['hidden']     = monochrome_logs ? '' : "\033[8m"

        // Regular Colors
        colorcodes['black']      = monochrome_logs ? '' : "\033[0;30m"
        colorcodes['red']        = monochrome_logs ? '' : "\033[0;31m"
        colorcodes['green']      = monochrome_logs ? '' : "\033[0;32m"
        colorcodes['yellow']     = monochrome_logs ? '' : "\033[0;33m"
        colorcodes['blue']       = monochrome_logs ? '' : "\033[0;34m"
        colorcodes['purple']     = monochrome_logs ? '' : "\033[0;35m"
        colorcodes['cyan']       = monochrome_logs ? '' : "\033[0;36m"
        colorcodes['white']      = monochrome_logs ? '' : "\033[0;37m"

        // Bold
        colorcodes['bblack']     = monochrome_logs ? '' : "\033[1;30m"
        colorcodes['bred']       = monochrome_logs ? '' : "\033[1;31m"
        colorcodes['bgreen']     = monochrome_logs ? '' : "\033[1;32m"
        colorcodes['byellow']    = monochrome_logs ? '' : "\033[1;33m"
        colorcodes['bblue']      = monochrome_logs ? '' : "\033[1;34m"
        colorcodes['bpurple']    = monochrome_logs ? '' : "\033[1;35m"
        colorcodes['bcyan']      = monochrome_logs ? '' : "\033[1;36m"
        colorcodes['bwhite']     = monochrome_logs ? '' : "\033[1;37m"

        // Underline
        colorcodes['ublack']     = monochrome_logs ? '' : "\033[4;30m"
        colorcodes['ured']       = monochrome_logs ? '' : "\033[4;31m"
        colorcodes['ugreen']     = monochrome_logs ? '' : "\033[4;32m"
        colorcodes['uyellow']    = monochrome_logs ? '' : "\033[4;33m"
        colorcodes['ublue']      = monochrome_logs ? '' : "\033[4;34m"
        colorcodes['upurple']    = monochrome_logs ? '' : "\033[4;35m"
        colorcodes['ucyan']      = monochrome_logs ? '' : "\033[4;36m"
        colorcodes['uwhite']     = monochrome_logs ? '' : "\033[4;37m"

        // High Intensity
        colorcodes['iblack']     = monochrome_logs ? '' : "\033[0;90m"
        colorcodes['ired']       = monochrome_logs ? '' : "\033[0;91m"
        colorcodes['igreen']     = monochrome_logs ? '' : "\033[0;92m"
        colorcodes['iyellow']    = monochrome_logs ? '' : "\033[0;93m"
        colorcodes['iblue']      = monochrome_logs ? '' : "\033[0;94m"
        colorcodes['ipurple']    = monochrome_logs ? '' : "\033[0;95m"
        colorcodes['icyan']      = monochrome_logs ? '' : "\033[0;96m"
        colorcodes['iwhite']     = monochrome_logs ? '' : "\033[0;97m"

        // Bold High Intensity
        colorcodes['biblack']    = monochrome_logs ? '' : "\033[1;90m"
        colorcodes['bired']      = monochrome_logs ? '' : "\033[1;91m"
        colorcodes['bigreen']    = monochrome_logs ? '' : "\033[1;92m"
        colorcodes['biyellow']   = monochrome_logs ? '' : "\033[1;93m"
        colorcodes['biblue']     = monochrome_logs ? '' : "\033[1;94m"
        colorcodes['bipurple']   = monochrome_logs ? '' : "\033[1;95m"
        colorcodes['bicyan']     = monochrome_logs ? '' : "\033[1;96m"
        colorcodes['biwhite']    = monochrome_logs ? '' : "\033[1;97m"

        return colorcodes
    }

    //
    // Does what is says on the tin
    //
    public static String dashedLine(monochrome_logs) {
        Map colors = logColours(monochrome_logs)
        return "${colors.dim}--------------------------------------------------------------------------------${colors.reset}"
    }

    // epi2me-labs logo
    public static String logo(workflow, monochrome_logs) {
        Map colors = NfcoreTemplate.logColours(monochrome_logs)
        String workflow_name = workflow.manifest.name.split("/")[1]
        String workflow_version = version(workflow)
        String.format(
            """
            ${colors.igreen}||||||||||   ${colors.reset}${colors.dim}_____ ____ ___ ____  __  __ _____
            ${colors.igreen}||||||||||  ${colors.reset}${colors.dim}| ____|  _ \\_ _|___ \\|  \\/  | ____|
            ${colors.yellow}|||||       ${colors.reset}${colors.dim}|  _| | |_) | |  __) | |\\/| |  _|
            ${colors.yellow}|||||       ${colors.reset}${colors.dim}| |___|  __/| | / __/| |  | | |__
            ${colors.iblue}||||||||||  ${colors.reset}${colors.dim}|_____|_|  |___|_____|_|  |_|_____|
            ${colors.iblue}||||||||||  ${colors.reset}${colors.bold}${workflow_name} ${workflow_version}${colors.reset}
            ${NfcoreTemplate.dashedLine(monochrome_logs)}
            """.stripIndent()
        )
    }
}


