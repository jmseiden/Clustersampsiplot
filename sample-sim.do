
cap program drop clustersampsiplot

program clustersampsiplot
		version 16
		syntax, mdes(real) rho(real) [clusters(numlist missingok >0)] [ns_per_cluster(numlist missingok >0)] [base_correl(real 0)] [alpha(real .05)] [beta(real .8)] [savesims(string)] [mdes2(numlist missingok max=1)] [rho2(numlist missingok max=1)] [base_correl2(numlist missingok max=1)] [alpha2(numlist missingok max=1)] [beta2(numlist missingok max=1)]
		
***Check to make sure arguments are OK and whether it is a single or double 

capture: confirm existence `clusters' 

	if _rc == 0 {
		
		capture: confirm existence `ns_per_cluster' 
		
			if _rc == 0 {

			display as error "Only clusters or ns_per_cluster may be specified"
			
			exit
			
			}

	}

capture: confirm existence `clusters' 

	if _rc != 0 {

		capture: confirm existence `ns_per_cluster' 
	
		display as error "Clusters or ns_per_cluster must be specified"
		
		exit
					
			}

***Check if the simulation will be with one or two sets of assumptions and set undefined secondary assumptions to original

	if "`mdes2'" != "" | "`rho2'" != "" | "`base_correl2'" != "" | "`alpha2'" != "" | "`beta2'" != "" {

		local double = 1

		if "`mdes2'" != "" {
			local mdes2 = `mdes2'
			} 
			else {
				local mdes2 = `mdes'
			}

		if "`rho2'" != "" {
			local rho2 = `rho2'
			} 
			else {
				local rho2 = `rho'
			}

		if "`base_correl2'" != "" {
			local base_correl2 = `base_correl2'
			} 
			else { 
				local base_correl2 = `base_correl'
			}
		if "`alpha2'" != "" {
			local alpha2 = `alpha2'
			} 
			else {
				local alpha2 = `alpha'
			}

		if "`beta2'" != "" {
			local beta2 = `beta2'
			} 
			else {
				local beta2 = `beta'
			}				
	}

	else {

		local double = 0 
	}

	*Calculate required ns_per_cluster for a given number of clusters
	capture: confirm existence `clusters'

	if _rc == 0 {

		display as result "Calculating ns_per_cluster for clusters: `clusters'"

		preserve
			
		clear all 
		
		*Single set of assumptions		
		if `double' == 0 {

			*Check if there might be a situation with no solution and report an error if so
			
			local min : word 1 of `clusters'

			forval j = 2/`: word count `clusters'' {
			    local min = min(`min', `: word `j' of `clusters'')
			}

			cap clustersampsi, mu1(`mdes') mu2(0) base_correl(`base_correl') alpha(`alpha') beta(`beta') rho(`rho') k(`min')

			local toosmall = r(min_k)

			if "`toosmall'" != "." {

				display as error "One or more cluster solutions could not be solved for. Try setting the minimum number of clusters above `min'"
			
			}

			*Extract the total number of elements and create a matrix of the correct dimensions

			local length: word count `clusters'

			matrix define results = J(`length',3,.)

			matrix colnames results = clusters npercluster totalsamplesize
			
			*Run the clustersampsi command to get the results and put them in the matrix
			
			local row = 1
			
			foreach x in `clusters' {
			
				qui: cap clustersampsi, mu1(`mdes') mu2(0) base_correl(`base_correl') alpha(`alpha') beta(`beta') rho(`rho') k(`x')

				local ns_per_cluster = r(m)

				matrix results[`row',1] = `x'

				matrix results[`row',2] = `ns_per_cluster'
				
				matrix results[`row',3] = `ns_per_cluster' * `x'
				
				local ++row

			}
			
			*Create the results and plot

			svmat results, names(col)

			la var npercluster "Required sample per cluster"
			la var clusters "Numbers of clusters (per arm)"
			la var totalsamplesize "Total sample size required"

			graph twoway (connected npercluster cluster), ///
					note("MDES=`mdes'; ICC=`rho'; Baseline correlation=`rho'; Power=`beta'; Alpha = `alpha'")

			*Save results if option is enabled
			capture: confirm existence `savesims' 
			
			if _rc == 0 {
			
				save "`savesims'", replace
			
			}
		}

		*Two sets of assumptions		
		if `double' == 1 {

			*Extract the total number of elements and create a matrix of the correct dimensions

			local length: word count `clusters'

			matrix define results = J(`length',5,.)

			matrix colnames results = clusters npercluster1 totalsamplesize1 npercluster2 totalsamplesize2
			
			*Run the clustersampsi command to get the results and put them in the matrix
			
			local row = 1
			
			foreach x in `clusters' {

				matrix results[`row',1] = `x'
			
				*First set of assumptions
				qui: cap clustersampsi, mu1(`mdes') mu2(0) base_correl(`base_correl') alpha(`alpha') beta(`beta') rho(`rho') k(`x')

				local ns_per_cluster1 = r(m)

				matrix results[`row',2] = `ns_per_cluster1'

				matrix results[`row',3] = `ns_per_cluster1' * `x'

				*Second set of assumptions
				di "mu1(`mdes2') mu2(0) base_correl(`base_correl2') alpha(`alpha2') beta(`beta2') rho(`rho2') k(`x')"
				qui: cap clustersampsi, mu1(`mdes2') mu2(0) base_correl(`base_correl2') alpha(`alpha2') beta(`beta2') rho(`rho2') k(`x')

				local ns_per_cluster2 = r(m)

				matrix results[`row',4] = `ns_per_cluster2'
				
				matrix results[`row',5] = `ns_per_cluster2' * `x'
				
				local ++row

			}
			
			*Create the results and plot

			svmat results, names(col)

			la var clusters "Numbers of clusters (per arm)"
			la var npercluster1 "Sample per cluster (MDES=`mdes'; ICC=`rho'; Base. corr.=`base_correl'; Power=`beta'; Alpha = `alpha')"
			la var totalsamplesize1 "Total sample size (MDES=`mdes'; ICC=`rho'; Base. corr.=`base_correl'=`rho'; Power=`beta'; Alpha = `alpha')"
			la var npercluster2 "Sample per cluster (MDES=`mdes2'; ICC=`rho2';  Base. corr.=`base_correl2'; Power=`beta2'; Alpha = `alpha2')"
			la var totalsamplesize2 "Total sample size (MDES=`mdes2'; ICC=`rho2';  Base. corr.=`base_correl2'; Power=`beta2'; Alpha = `alpha2')"

			graph twoway (connected npercluster1 clusters) (connected npercluster2 cluster), ///
				ytitle("Required sample size per cluster") legend(ring(0))

			*Save results if option is enabled
			capture: confirm existence `savesims' 
			
			if _rc == 0 {
			
				save "`savesims'", replace
			
			}
		}		
				
		restore		
		
		exit

	}

	*Calculate required ns_per_cluster for a given number of clusters

	capture: confirm existence `ns_per_cluster' 

		if _rc == 0 {
		
			display as result "Calculating clusters for ns_per_cluster: `ns_per_cluster'"
			
			*Calculate required clusters for a given number of ns_per_cluster

			preserve
			
			clear all 

			*Code for single simulation
			
			if `double' == 0 {			

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
				
			*Two sets of assumptions		
			if `double' == 1 {

				*Extract the total number of elements and create a matrix of the correct dimensions

				local length: word count `ns_per_cluster'

				matrix define results = J(`length',5,.)

				matrix colnames results = nperclusters clusters1 totalsamplesize1 clusters2 totalsamplesize2
				
				*Run the clustersampsi command to get the results and put them in the matrix
				
				local row = 1
				
				foreach x in `ns_per_cluster' {

					matrix results[`row',1] = `x'
				
					*First set of assumptions
					qui: cap clustersampsi, mu1(`mdes') mu2(0) base_correl(`base_correl') alpha(`alpha') beta(`beta') rho(`rho') m(`x')

					local clusters1 = r(m)

					matrix results[`row',2] = `clusters1'

					matrix results[`row',3] = `clusters1' * `x'

					*Second set of assumptions
					di "mu1(`mdes2') mu2(0) base_correl(`base_correl2') alpha(`alpha2') beta(`beta2') rho(`rho2') k(`x')"
					qui: cap clustersampsi, mu1(`mdes2') mu2(0) base_correl(`base_correl2') alpha(`alpha2') beta(`beta2') rho(`rho2') m(`x')

					local clusters2 = r(m)

					matrix results[`row',4] = `clusters2'
					
					matrix results[`row',5] = `clusters2' * `x'
					
					local ++row

				}
				
				*Create the results and plot

				svmat results, names(col)

				la var clusters "Numbers of clusters (per arm)"
				la var npercluster1 "Sample per cluster (MDES=`mdes'; ICC=`rho'; Base. corr.=`base_correl'; Power=`beta'; Alpha = `alpha')"
				la var totalsamplesize1 "Total sample size (MDES=`mdes'; ICC=`rho'; Base. corr.=`base_correl'=`rho'; Power=`beta'; Alpha = `alpha')"
				la var npercluster2 "Sample per cluster (MDES=`mdes2'; ICC=`rho2';  Base. corr.=`base_correl2'; Power=`beta2'; Alpha = `alpha2')"
				la var totalsamplesize2 "Total sample size (MDES=`mdes2'; ICC=`rho2';  Base. corr.=`base_correl2'; Power=`beta2'; Alpha = `alpha2')"

				graph twoway (connected npercluster1 clusters) (connected npercluster2 cluster), ///
					ytitle("Required sample size per cluster") legend(ring(0))

				*Save results if option is enabled
				capture: confirm existence `savesims' 
				
				if _rc == 0 {
				
					save "`savesims'", replace
				
				}
			}		

				restore

				exit
						
				}
			}


end

