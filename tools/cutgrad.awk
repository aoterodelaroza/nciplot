#!/usr/bin/awk -f
##############################################################################
# A script to discard data with a reduced density gradient less than         #
# gradmin or greater than gradmax                                            #
#                                                                            #
# Usage:                                                                     #
#  user:$ ./cutgrad.awk gradmin gradmax file.dat                             #
#                                                                            #
##############################################################################

BEGIN{
    if (ARGC < 3){
	print "Usage: cutgrad.awk gradmin gradmax file.dat"
	exit
    }
    #Read reduced density gradient range
    gradmin = ARGV[1]; gradmax = ARGV[2]
    #Remove range from arguments list
    for (i=1;i<=ARGC-1;i++)
	ARGV[i] = ARGV[i+2]
    ARGC -= 2
}
#Print data in the range (gradmin < grad < gradmax) to the screen
($2<gradmax)&&($2>gradmin)
