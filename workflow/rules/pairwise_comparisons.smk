pw_dirs = []
if pairs:
	dirs = ['increased', 'decreased', 'unchanged']
	for i in range(3):
		for j in range(3):
			pw_dirs.extend([dirs[i] + "_" + dirs[j]])


rule prepare_pairwise_results:
	input:
		blacklist = os.path.join(annotation_path, "blacklist.rds"),
        features = os.path.join(annotation_path, "features.rds"),
		greylist =  os.path.join(grey_path, "greylists.rds"),
        gtf_gene = os.path.join(annotation_path, "gtf_gene.rds"),
        hic = os.path.join(annotation_path, "hic.rds"),
		regions = os.path.join(annotation_path,"gene_regions.rds"),
		results1 = os.path.join(
			diff_path, "{tgt1}", "{tgt1}_{comp1}-differential-signal.rds"
		),
		results1 = os.path.join(
			diff_path, "{tgt2}", "{tgt2}_{comp2}-differential-signal.rds"
		),
		script = os.path.jon("workflow", "scripts", "pairwise_reuslts.R"),
		yaml = "config/params.yml"
	output:
		rds = os.path.join(
			pairs_path, "{tgt1}_{comp1}_{tgt2}_{comp2}", 
			"{tgt1}_{comp1}_{tgt2}_{comp2}-pairwise_results.rds"
		),
		bed = expand(
			os.path.join(
				pairs_path, "{{tgt1}}_{{comp1}}_{{tgt2}}_{{comp2}}",
				"{{tgt1}}_{{comp1}}_{{tgt2}}_{{comp2}}-{f}.bed.gz"
			),
			f = pw_dirs
		)
	params:
		config['pairwise']['default'] # Change later
	threads: 2
	conda: "../envs/rmarkdown.yml"
	resources:
		runtime = "15m"
		mem_mb = 16000
	log: os.path.join(log_path, "prepare_pairwise_results", "{tgt1}_{comp1}_{tgt2}_{comp2}.log")
	script:
        "../scripts/pairwise_results.R"		

# rule create_pairwise_comparisons_rmd:
# 	input:
# 		chk = ALL_CHECKS,
# 		module_pw = "workflow/modules/pairwise_comparison.Rmd",
# 		r = "workflow/scripts/create_pairwise_comparison.R"
# 	output:
# 		rmd = expand(
# 			os.path.join(
# 				"analysis", 
# 				"{{t1}}_{{ref1}}_{{treat1}}_{{t2}}_{{ref2}}_{{treat2}}_{f}"
# 				),
# 			f = "pairwise_comparison.Rmd"
# 		)
# 	params:
# 		threads = 4,
# 	conda: "../envs/rmarkdown.yml"
# 	threads: 1
# 	resources:
# 		runtime = "1m",
# 		mem_mb = 512
# 	log: log_path + "/create_rmd/create_{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}_pairwise_comparison_rmd"
# 	shell:
# 		"""
# 		## Create the generic markdown header
#         Rscript --vanilla \
#             {input.r} \
#             {wildcards.t1} \
#             {wildcards.ref1} \
#             {wildcards.treat1} \
# 			{wildcards.t2} \
#             {wildcards.ref2} \
#             {wildcards.treat2} \
#             {params.threads} \
#             {output.rmd} &>> {log}

#         ## Add the remainder of the module as literal text
#         cat {input.module_pw} >> {output.rmd}
# 		"""

# rule compile_pairwise_comparisons_html:
# 	input:
# 		annotations = ANNOTATION_RDS,
# 		blacklist = blacklist,
# 		config = "config/config.yml",
# 		module_rna = "workflow/modules/rnaseq_pairwise.Rmd",
# 		results_t1 = os.path.join(
# 			"docs", "{t1}_{ref1}_{treat1}_differential_signal.html"
# 		),
# 		results_t2 = os.path.join(
# 			"docs", "{t2}_{ref2}_{treat2}_differential_signal.html"
# 		),
# 		rmd = expand(
# 			os.path.join(
# 				"analysis", 
# 				"{{t1}}_{{ref1}}_{{treat1}}_{{t2}}_{{ref2}}_{{treat2}}_{f}"
# 				),
# 			f = "pairwise_comparison.Rmd"
# 		),
# 		rmd_config = "config/rmarkdown.yml",
# 		scripts = os.path.join("workflow", "scripts", "custom_functions.R"),
# 		setup = rules.create_setup_chunk.output,
# 		yaml = rules.create_site_yaml.output
# 	output:
# 		html = expand(
# 			os.path.join(
# 				"docs",
# 				"{{t1}}_{{ref1}}_{{treat1}}_{{t2}}_{{ref2}}_{{treat2}}_{f}"
# 			),
# 			f = "pairwise_comparison.html"
# 		),
# 		fig_path = directory(
# 			os.path.join(
# 				"docs",
# 				"{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}_pairwise_comparison_files"
# 			)
# 		),
# 		csv = expand(
# 			os.path.join(
# 				"output", "pairwise_comparisons", "{{t1}}_{{t2}}",
# 			"{{t1}}_{{ref1}}_{{treat1}}-{{t2}}_{{ref2}}_{{treat2}}-{f}"
# 			),
# 			f = [
# 				'pairwise_comparison.csv.gz', 'enrichment.csv', 
# 				'rnaseq_enrichment.csv'
# 				]
# 		),
# 		rds = os.path.join(
# 			"output", "pairwise_comparisons", "{t1}_{t2}",
# 			"{t1}_{ref1}_{treat1}-{t2}_{ref2}_{treat2}-all_windows.rds"
# 		),
# 		renv = temp(
# 			os.path.join(
# 				"output", "envs",
# 				"{t1}_{ref1}_{treat1}-{t2}_{ref2}_{treat2}-pairwise_comparison.RData"
# 			)
# 		)
# 	params:
# 		asset_path = os.path.join(
# 			"docs", "assets", "{t1}_{ref1}_{treat1}-{t2}_{ref2}_{treat2}"
# 		)
# 	conda: "../envs/rmarkdown.yml"
# 	threads: 4
# 	resources:
# 		runtime = "2h",
# 		mem_mb = 8192
# 	log: "workflow/logs/pairwise/{t1}_{ref1}_{treat1}_{t2}_{ref2}_{treat2}_pairwise_comparison.log"
# 	resources:
# 		mem_mb = 32768,
# 		runtime = "2h"
# 	shell:
# 		"""
#         R -e "rmarkdown::render_site('{input.rmd}')" &>> {log}
# 		"""