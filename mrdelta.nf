#!/usr/bin/env nextflow


/*
MRdelta all vs all

score:  [tc, sp, col, sspa, pastml2, pastml3]
metric: [homo, whomo, whomo2, len, ngap, ngap2, seqid, pastml]  - negative
metric: [tc, sp, col, sspa, pastml2, pastml3]  - positive


score:  [homo, whomo, whomo2, len, ngap, ngap2, seqid, pastml]
metric: [tc, sp, col, sspa, pastml2, pastml3]  - negative
metric: [homo, whomo, whomo2, len, ngap, ngap2, seqid, pastml]  - positive
*/

params.score = ["sp"]
// params.metric = ["homo", "ngap", "pastml1000", "wpastml1000", "hpastml1000", "fragmentation"]
// params.metric = ["tc"] 
params.metric = ["fragmentation_avg", "fragmentation_sum"]


params.mrdelta = [-1, 0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0]
params.maintain = "bucket aligner family"
params.deltaby = "tree"
params.norm = "" // "-norm avgLenXnseq"
params.corr = "-corr negative" 
params.maxseq = "" // "-maxseq 1000"


params.csv = "/users/cn/sjin/projects/homoplasy/pastml/pastml_parent_1000_CLUSTALO/homoplasy_pastml_fragmentation.csv"
params.output = "/users/cn/sjin/projects/homoplasy/pastml/pastml_parent_1000_CLUSTALO/mrdelta_fragmentation_below1000"



log.info """\
         Running MRdelta - python version"
         ======================================="
         Input csv                                          : ${params.csv}
         Evaluation:
            Scores                                          : ${params.score}
            Metrics                                         : ${params.metric}
            Mrdelta                                         : ${params.mrdelta}
            Between                                         : ${params.deltaby}
         Output directory                                   : ${params.output}

         """
         .stripIndent()




process mrdelta {
    
    cache false

    tag "${score}_${metric}_mrdelta${mrdelta}"
    publishDir "${params.output}/mrdelta${mrdelta}", mode: 'copy', overwrite: true

    input:
        path(csv) from params.csv
        each score from params.score
        each metric from params.metric
        each mrdelta from params.mrdelta
        val(maintain) from params.maintain
        val(deltaby) from params.deltaby
        val(norm) from params.norm
        val(corr) from params.corr
        val(maxseq) from params.maxseq
    
    output:
        file("${score}_${metric}_mrdelta${mrdelta}*.tsv") into mrdeltaOut

    when:
        if(score == metric){false}else{true}

    script:
        """
        python ${baseDir}/bin/mrdelta.py -csv ${csv} \
            -score ${score} -metric ${metric} -mrdelta ${mrdelta} -outdir . -maintain ${maintain} -deltaby ${deltaby} ${norm} ${corr} ${maxseq}

        """
}