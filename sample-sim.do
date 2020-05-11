
cap program drop clustersampsiplot

program clustersampsiplot
		version 16
		syntax, mdes(real) rho(real) [clusters(numlist missingok >0)] [ns_per_cluster(numlist missingok >0)] [base_correl(real max = 1 min = 0 0)] [alpha(real .05)] [beta(real .8)] [savesims(string)] mdes2(numlist max = 1)
		
		
***Check to make sure arguments are OK

capture: confirm existence `clusters' 

	if _rc == 0 {
		
		capture: confirm existence `ns_per_cluster' 
		
			if _rc == 0 {

			display as error "Only clusters or ns_per_cluster may be specified"
			
			exit
			
			}
		
		display as result "Calculating ns_per_cluster for clusters: `clusters'"
		
***Check if the simulation will be with one or two sets of assumptions

	
	*Calculate required ns_per_cluster for a given number of clusters
	
		preserve
		
		clear all 

		*Extract the total number of elements and create a matrix of the correct dimensions

		local length: word count `clusters'

		matrix define results = J(`length',3,.)

		matrix colnames results = npercluster cluster totalsamplesize
		
		*Run the clustersampsi command to get the results and put them in the matrix
		
		local row = 1
		
		foreach x in `clusters' {
		
			qui: cap clustersampsi, mu1(`mdes') mu2(0) base_correl(`base_correl') alpha(`alpha') beta(`beta') rho(`rho') k(`x')

			local ns_per_cluster = r(m)

			matrix results[`row',1] = `ns_per_cluster'

			matrix results[`row',2] = `x'
			
			matrix results[`row',3] = `ns_per_cluster' * `x'
			
			local ++row

		}
		
		*Create the results and plot

		svmat results, names(col)

		la var npercluster "Required sample per cluster"
		la var cluster "Numbers of clusters (per arm)"
		la var totalsamplesize "Total sample size required"

		graph twoway (connected npercluster cluster), ///
				note("Minimum Detectable Effect Size: `mdes'; Intra-class correlation: `rho'; Power: `beta'; Error rate: `alpha'")

		*Save results if option is enabled
		capture: confirm existence `savesims' 
		
		if _rc == 0 {
		
			save "`savesims'", replace
		
		}
				
		restore		
		
		exit

	}


capture: confirm existence `ns_per_cluster' 

	if _rc == 0 {
	
		display as result "Calculating clusters for ns_per_cluster: `ns_per_cluster'"
		
		*Calculate required clusters for a given number of ns_per_cluster

		preserve
		
		clear all 

			*Extract the total number of elements and create a matrix of the correct dimensions

			local length: word count `ns_per_cluster'

			matrix define results = J(`length',3,.)

			matrix colnames results = npercluster cluster totalsamplesize
			
			*Run the clustersampsi command to get the results and put them in the matrix
			
			local row = 1

			foreach x in `ns_per_cluster' {

				qui: cap clustersampsi, mu1(`mdes') mu2(0) base_correl(`base_correl') alpha(`alpha') beta(`beta') rho(`rho') m(`x')

				local cluster = r(k)

				matrix results[`row',1] = `x'

				matrix results[`row',2] = `cluster'
				
				matrix results[`row',3] = `cluster' * `x'

				local ++row

			}
			
			*Create the results and plot

			svmat results, names(col)
			
			la var npercluster "Sample size per cluster"
			la var cluster "Required numbers of clusters (per arm)"
			la var totalsamplesize "Total sample size required"

			graph twoway (connected cluster npercluster), ///
				note("Minimum Detectable Effect Size: `mdes'; Intra-class correlation: `rho'; Power: `beta'; Error rate: `alpha'")
			
			*Save results if option is enabled
			capture: confirm existence `savesims' 
			
			if _rc == 0 {
			
				save "`savesims'", replace
			
			}	
			
			exit
					
			}
	
	if _rc != 0 {
	
		display as error "Clusters or ns_per_cluster must be specified"
		
		exit
					
			}

	}

end

