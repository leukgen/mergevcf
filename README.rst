Merge SV VCFs
=============

This set of routines merges SV VCF files by calls, labeling calls by the callers that
made them in an INFO field, Callers=.   Most of the work goes into normalizing
SV calls, which are treated as paired of breakpoints that are equal if they fall within
a user-specified window.

------------
Known Issues
------------

* INFO and FORMAT merged for multiple from same caller with only be for last record
* Individual Caller Genotype not kept
