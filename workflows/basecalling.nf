include {
    create_metamap;
} from '../lib/bamingress'

include { wf_dorado } from './_basecalling'

def create_bam_channel(def merge_out){
    merge_out.map{
        // map basecalled cram to (xam_path, xam_index, xam_meta) tuple
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
        // returns: (cram), (crai), map to ([cram, crai])
        out = wf_dorado(input_path, ref)
        bam_channel = create_bam_channel(out.cram.concat(out.crai).toList()).view()
    emit:
        bam_channel
}
