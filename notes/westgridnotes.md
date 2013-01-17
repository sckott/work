# westgrid notes

## notes
+ avoid creating many small files
	+ what is the quota?
		+ 
		+ 
+ system uses Lustre file system
+ TORQUE and Moab commands
	+ qsub, -l procs=48 will distributed jobs across 48 nodes
	+ showq
	+ showstart jobid - show me the estimated start time of the jobid
	+ showbf - show the back-fill window
	+ qstat
	+ qdel
+ Can apply for a larger resource allocation [here](http://www.westgrid.ca/support/accounts/resource_allocations)
+ can ask for walltime extensions to westgrid support

## Questions
+ I need example job scripts for R and python jobs, in parallel and not, etc.
+ 	+ they didn't really give any
+ Can I run web based RStudio from westgrid?
	+ no, for now - probably not ever
	+ they have done it for some web interfaces, but has to be investigated
+ R
	+ different versions of R - 
		+ can install on user directory, and won't affect global system
	+ they are reluctant to upgrade to new versions or R
	+ Can try to use MPI to parallelize R scripts