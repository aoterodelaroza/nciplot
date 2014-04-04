#! /usr/bin/awk -f

BEGIN{
    if (ARGC < 3){
	print "Usage: cutden.awk rhomin rhomax [file1 file2 ...]"
	exit
    }
    rhomin = ARGV[1]; rhomax = ARGV[2]
    for (i=1;i<=ARGC-1;i++)
	ARGV[i] = ARGV[i+2]
    ARGC -= 2
}
($1<rhomax)&&($1>rhomin)
