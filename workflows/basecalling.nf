include {
    create_metamap;
} from '../lib/bamingress'

include { wf_dorado } from './_basecalling'

def create_bam_channel(def merge_out){
    // map basecalled cram to (xam_path, xam_index, xam_meta) tuple
    merge_out.map{
        it -> tuple(it[0], it[1], create_metamap([
            output: true, // write this to out_dir
            is_cram: true, // we know we basecalled it to cram
        ]))
    }
}

// shim the wf-basecalling workflow to emit the "bam channel" as if this had
// passed through humvars lib/bamingress
workflow basecalling {
    take:
        input_path
        ref
    main:
        // pass args straight to actual basecalling workflow
        crams = wf_dorado(
            input_path,
            ref,
            params.basecaller_cfg, params.basecaller_model_path,
            params.remora_cfg, params.remora_model_path,
        )
        // annotate results and emit
        pass = create_bam_channel(crams.pass)
        fail = create_bam_channel(crams.fail)
    emit:
        pass = pass
        fail = fail
}
