.. package_inventory:
.. _package_inventory:

****************
Package Manifest
****************

The `GitHub repository <https://github.com/>`_ checkout should contain the 
following files:

root
	:ref:`Annotations <annotations>`

		:ref:`ERCC_Controls_Analysis.txt <ERCC_Controls_Analysis_tsvdoc>`
				
		:ref:`ERCC_Controls_Annotation.txt <ERCC_Controls_Annotation_tsvdoc>`
		
		:ref:`ERCC92.gtf <ERCC92_tsvdoc>`
		
		:ref:`get_gene_info.pl <get_gene_info_perldoc>`

		:ref:`RNA_sample_details.tsv <RNA_sample_details_tsvdoc>`
		
		:ref:`Saccharomyces_cerevisiae_with_spike_ins.EF4.68.gtf <Saccharomyces_cerevisiae_with_spike_ins_tsvdoc>`
		
		:ref:`Saccharomyces_cerevisiae.EF4.68.gtf <Saccharomyces_cerevisiae_tsvdoc>`
	
	
	:ref:`Bad_replicate_identification <badreps>`
		
		:ref:`bad_replicates.pl <bad_replicates_perldoc>`
		
		:ref:`compare_replicates.pl <compare_replicates_perldoc>`
	
	
	:ref:`Bootstrapping <bootstrapping>`
		
		:ref:`add_gene_name_column.pl <add_gene_name_column_perldoc>`
		
		:ref:`combine_bs_pvalues.pl <combine_bs_pvalues_perldoc>`
		
		:ref:`exclude_badreps.lst <exclude_badreps_tsvdoc>`
		
		:ref:`generic_wrapper.py <generic_wrapper_autodoc>`
		
		:ref:`grid_launcher_DE.pl <grid_launcher_DE_perldoc>`
		
		:ref:`grid_launcher_powertest_cuffdiff.pl <grid_launcher_powertest_cuffdiff_perldoc>`
				
		:ref:`grid_launcher_powertest.pl <grid_launcher_powertest_perldoc>`
		
		:ref:`group_by_gene.pl <group_by_gene_perldoc>`
		
		:ref:`make_cuffdiff_db.pl <make_cuffdiff_db_perldoc>`
		
		:ref:`one_bootstrap_ET93.pl <one_bootstrap_ET93_perldoc>`
		
		:ref:`one_bs_pvalue.pl <one_bs_pvalue_perldoc>`
		
		:ref:`one_permutest.pl <one_permutest_perldoc>`
		
		:ref:`one_powertest_simple.pl <one_powertest_simple_perldoc>`

		:ref:`qsub_bootstraps_withexcludelist.sh <qsub_bootstraps_withexcludelist_shdoc>`
		
		:ref:`qsub_bootstraps_withoutexludelist.sh <qsub_bootstraps_withoutexludelist_shdoc>`
		
		:ref:`qsub_WTfdr_bootstraps.sh <qsub_WTfdr_bootstraps_shdoc>`

		:ref:`simple_de_test.pl <simple_de_test_perldoc>`
	
			
	:ref:`DE_tool_comparison <DEtoolcomp>`
	
		:ref:`gene_ranking_sigprop.pl <gene_ranking_sigprop_perldoc>`
		
		:ref:`gene_ranking.pl <gene_ranking_perldoc>`
		
		:ref:`make_powerstats_db.pl <make_powerstats_db_perldoc>`
		
		:ref:`one_bs_powerstats.pl <one_bs_powerstats_perldoc>`
	
	
	:ref:`DE_tool_scripts <DEtools>`

		:ref:`bayseq.R <bayseq_Rdoc>`

		:ref:`cuffdiff_pq.pl <cuffdiff_pq_perldoc>`
		
		:ref:`degseq.R <degseq_Rdoc>`

		:ref:`deseq.R <deseq_Rdoc>`
		
		:ref:`deseq2.R <deseq2_Rdoc>`
		
		:ref:`ebseq.R <ebseq_Rdoc>`
		
		:ref:`edgeR.R <edgeR_Rdoc>`
				
		:ref:`limma.R <limma_Rdoc>`
		
		:ref:`noiseq.R <noiseq_Rdoc>`
		
		:ref:`poissonseq.R <poissonseq_Rdoc>`
		
		:ref:`run_cuffdiff.pl <run_cuffdiff_perldoc>`
		
		:ref:`samseq.R <samseq_Rdoc>`

		:ref:`t-test_equal.R <ttest_equal_Rdoc>`
		
		:ref:`t-test_unequal.R <ttest_unequal_Rdoc>`
	
	
	:ref:`Distribution_tests <Dists>`

		:ref:`distribution_test.pl <distribution_test_perldoc>`
		
		:ref:`grid_launcher_nbtest.pl <grid_launcher_nbtest_perldoc>`
		
		:ref:`one_nb_test.pl <one_nb_test_perldoc>`
		
		:ref:`poisson_reps.pl <poisson_reps_perldoc>`
		
	:ref:`General scripts and files <General>`

		:ref:`combine_replicates.pl <combine_replicates_perldoc>`
		
		:ref:`de_tests.txt <de_tests_configfile>`
		
		:ref:`de_tests_samecond.txt <de_tests_samecond_configfile>`
		
		:ref:`defs.dat <defs_dat_configfile>`
		
		:ref:`exclude_badreps.lst <exclude_badreps_lst_configfile>`
		
		:ref:`genlist.tsv <genlist_configfile>`
		
		
	:ref:`Mapping <mapping>`
		
		:ref:`build_pileup.pl <build_pileup_perldoc>`
		
		:ref:`concat_fastq.pl <concat_fastq_perldoc>`
				
		:ref:`count_fastq_reads_lanes.pl <count_fastq_reads_lanes_perldoc>`
		
		:ref:`count_fastq_reads.pl <count_fastq_reads_perldoc>`
				
		:ref:`make_bamlinks.pl <make_bamlinks_perldoc>`
		
		:ref:`pileup_var.pl <pileup_var_perldoc>`
		
		:ref:`run_alignments.pl <run_alignments_perldoc>`
				
		:ref:`run_biological_alignments.pl <run_biological_alignments_perldoc>`
		
		:ref:`sort_and_index_bams.py <sort_and_index_bams_autodocs>`
		
		:ref:`summary_stats.pl <summary_stats_perldoc>`

		:ref:`unique_fastq_reads.pl <unique_fastq_reads_perldoc>`
		
	
	:ref:`Modules <modules>`

		:ref:`BootstrapDB.pm <BootstrapDB_perlmod>`
		
		:ref:`Cluster.pm <Tools_perlmod>`
		
		:ref:`Distribution.pm <Distribution_perlmod>`
		
		:ref:`GetOpt.pm <GetOpt_perlmod>`
		
		:ref:`GRNASeq.pm <GRNASeq_perlmod>`
		
		:ref:`housekeeping.py <housekeeping_autodoc>`
		
		:ref:`Powertests.pm <Powertests_perlmod>`
		
		:ref:`Tools.pm <Tools_perlmod>`
		
		:ref:`TreeView.pm <TreeView_perlmod>`

	
	:ref:`Normalization <normalisation>`

		:ref:`levelling_lanes.pl <levelling_lanes_perldoc>`
		
		:ref:`levelling.pl <levelling_perldoc>`
		
		:ref:`norms.pl <norms_perldoc>`
		
		:ref:`one_levelling.pl <one_levelling_perldoc>`

	
	:ref:`Plotting <plotting>`

		:ref:`cluster_replicates_pv.R <cluster_replicates_pv_Rdoc>`
		
		:ref:`cluster_tools_pv.R <cluster_tools_pv_Rdoc>`
		
		:ref:`compare_null.pl <compare_null_perldoc>`
		
		:ref:`compare_powerstats_db.pl <compare_powerstats_db_perldoc>`
		
		:ref:`plot_distest.pl <plot_distest_perldoc>`
		
		:ref:`plot_fcp.pl <plot_fcp_perldoc>`
				
		:ref:`plot_gene_examples.pl <plot_gene_examples_perldoc>`
		
		:ref:`plot_gene_pileup.pl <plot_gene_pileup_perldoc>`
				
		:ref:`plot_ms.pl <plot_ms_perldoc>`
				
		:ref:`plot_powerstats_db.pl <plot_powerstats_db_perldoc>`
		
		:ref:`plot_spikein_groups.pl <plot_spikein_groups_perldoc>`


	:ref:`Pre-processed Data <preprocessed_data>`
	
		:ref:`Gold-=standard Results <gold_standard_results>`
		
		:ref:`ENA_sample_mapping.pl <ENA_sample_mapping_perldoc>`
		
		:ref:`ENAdata_ERP004763_sample_mapping.tsv <ENAdata_ERP004763_sample_mapping_tsv>`
		
		:ref:`WT_countdata.tar.gz & WT_countdata.tar.gz <countdata.tar.gz>`

..	toctree::
	:maxdepth: 1
	