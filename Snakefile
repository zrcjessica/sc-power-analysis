configfile: "config.yml"

rule all:
    input: expand(config['out']+"/plots/{sample}_power_curves.pdf", sample = config['sample'])

rule est_params:
    input:
        counts = config['data'][config['sample']]['counts']
    output: config['out'] + "/{sample}_nb_params.csv"
    conda: "environment.yml"
    shell:
        "Rscript scripts/estimate_nb_params.R --count_mtx {input.counts} --out {output}"

rule power_analysis:
    input:
        counts = config['data'][config['sample']]['counts'],
        nb_params = rules.est_params.output
    params:
        alpha = config['params']['alpha'],
        nGenes = config['params']['nGenes'],
        pctDegs = config['params']['pctDegs'],
        test = config['params']['test'],
        nSims = config['params']['nSims']
    output: config['out'] + "/{sample}_{repls}repls_{cells}cells_{effect}effect.txt"
    conda: "environment.yml"
    shell:
        "Rscript scripts/power_analysis.R --count_mtx {input.counts} "
        "--params {input.nb_params} --out {output} "
        "--replicates {wildcards.repls} --effectSize {wildcards.effect} --nCells {wildcards.cells} "
        "--alpha {params.alpha} --nGenes {params.nGenes} " 
        "--pctDiff {params.pctDegs} --test {params.test} --nSims {params.nSims} "

rule aggregate_outputs:
    input: expand(rules.power_analysis.output[0], sample = '{sample}', \
        repls = config['params']['replicates'], \
        cells = config['params']['cells'], \
        effect = config['params']['effects'])
    output: config['out'] + "/{sample}_power_analysis_results.csv"
    conda: "environment.yml"
    shell: 
        "cat {input} > {output} &&"
        "rm {input}"

rule power_curves:
    input: rules.aggregate_outputs.output
    params: config['params']['sample_variable']
    output: config['out'] + "/plots/{sample}_power_curves.pdf"
    conda: "environment.yml"
    shell:
        "Rscript scripts/plot_power_curves.R --power_analysis {input} "
        " --sample_variable {params} --out {output}"
