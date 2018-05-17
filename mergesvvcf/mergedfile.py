from __future__ import print_function
import os
import pysam
import mergesvvcf.variantdict as variantdict

def mapped_to_chromosome(chrom):
    """
    Returns true if mapped to, eg, chr1 or X;
    false if mapped to other contig, eg GL*, MT*, hs*, M*
    """
    if chrom[0:2] in ["GL", "MT", "hs", "NC"] or chrom[0:1] == "M":
        return False
    return True

def int_if_possible(val):
    """
    Returns integer value of s if conversion succeeds, else s.
    """
    if type(val) is list:
        val = val[0]
    try:
        i = int(val)
    except ValueError:
        i = val
    return i

def make_info_dict(records, pos1, pos2):
    """ generate 'median' results for info fields from all records """
    fields = ['CHR2', 'END', 'SVTYPE', 'SVLEN']
    info = {}
    for field in fields:
        answers = []
        for record in records:
            if field in record.info:
                answers.append(int_if_possible(record.info[field]))
        nanswers = len(answers)
        if nanswers > 0:
            sorted_answers = sorted(answers)
            # doesn't quite deal with even #s correctly - can't average strings
            median_pos = int(nanswers / 2)
            median_answer = sorted_answers[median_pos]
            if not median_answer == 0:
                info[field] = median_answer

    if 'SVTYPE' in info and info['SVTYPE'] in ['DUP', 'DUP:TANDEM', 'DEL', 'INV']:
        if not 'SVLEN' in info:
            info['SVLEN'] = pos2-pos1
    return info

def bkptRefAltFromPair(loc1, loc2, refstr="N"):
    alt_after = loc1.__right__ == False

    if loc2 is None or loc2.__chrom__ is None:
        bkptstr = "."
    else:
        if loc2.__right__:
            delim = "["
        else:
            delim = "]"

        assert loc2.__strand__ == (loc1.__right__ != loc2.__right__)
        bkptstr = "%s%s:%d%s" % (delim, loc2.__chrom__, loc2.__pos__, delim)

    if alt_after:
        altstr = "%s%s" % (refstr, bkptstr)
    else:
        altstr = "%s%s" % (bkptstr, refstr)

    return refstr, altstr


def merge(filenames, programs, forceSV, outfile, slop=0, verbose=True,
        output_ncallers=False, min_num_callers=0,
        filterByChromosome=True, noFilter=False, debug=False):
    """Merge several VCFs from different programs into a new VCF file."""

    outvcf = pysam.VariantFile(os.devnull, 'w')
    if min_num_callers > 0:
        outvcf.header.filters.add("LOWSUPPORT", None, None, "Not called by enough callers in ensemble")

    # Returns true if the variant is PASS in the VCF file
    def passed_variant(record):
        """Did this variant pass?"""
        return "PASS" in list(record.filter) or len(list(record.filter)) == 0 or noFilter

    def infoString(callers, infodict):
        """
        Generate an INFO string from the INFO dictionary plus
        the list of callers.
        """
        info = {}
        fields = ['SVTYPE', 'SVLEN']
        for field in fields:
            if field in infodict:
                res = infodict[field]
                if type(res) is list:
                    res = res[0]
                info[field] = str(res)
        if output_ncallers:
            info["NumCallers"] = int(len(set(callers)))
        info["Callers"] = ",".join(list(set(callers)))
        return info

    calldict = variantdict.variantmap(awindow=0, svwindow=slop)
    for (infile, program) in zip(filenames, programs):
        count = 0
        try:
            vcf_reader = pysam.VariantFile(infile)

            for contig in sorted(vcf_reader.header.contigs):
                try:
                    outvcf.header.contigs.add(contig)
                except ValueError:
                    pass

            for record in vcf_reader.fetch():

                if not passed_variant(record):
                    continue

                if filterByChromosome and not mapped_to_chromosome(record.chrom):
                    continue

                if verbose:
                    if count == 0:
                        print(record, program)
                    count += 1
                    if count == 100:
                        count = 0

                calldict.addrecord(record, program, forceSV)
        except (RuntimeError, TypeError, NameError, AttributeError) as e:
            if debug:
                raise(e)
            else:
                pass

        outvcf.header.add_meta(
            key="caller",
            items=[
                ("ID", program),
                ("File", infile)
                ]
        )

    # Write the results in a master vcf file for the sample
    outvcf.header.info.add("Callers", ".", "String", "Callers that made this call")

    if output_ncallers:
        outvcf.header.info.add("NumCallers", 1, "Integer", "Number of callers that made this call")

    outvcf.header.info.add("SVTYPE", "1", "String", "Type of structural variant")
    outvcf.header.info.add("SVLEN", "1", "String", "Length of structural variant")

    print(outvcf.header, end="", file=outfile)

    for variant in sorted(calldict):
        callers = variant[2]
        num_callers = len(set(callers))
        passes = num_callers >= min_num_callers
        filterstring = "." if passes else "LOWSUPPORT"

        if len(variant) == 3:   # snv/indel
            loc, allele, callers = variant
            if allele is None:
                print("Allele is none: loc, allele, callers = ", loc, allele, callers)
                continue
            chrom, pos, _, _ = loc.asTuple()
            vcfrec = outvcf.new_record()
            vcfrec.contig = chrom
            vcfrec.start = pos
            vcfrec.ref = allele[0]
            vcfrec.alts = [allele[1]]
            vcfrec.filter.add(filterstring)
            vcfrec.info.__setitem__("Callers", ",".join(callers))
            if output_ncallers:
                vcfrec.info.__setitem__("NumCallers", str(len(set(callers))))
            print(vcfrec, end="", file=outfile)
        else:
            loc1, loc2, callers, medianPos1, medianPos2, recordscalled = variant
            if filterByChromosome and not mapped_to_chromosome(loc2.chrom):
                continue

            records = [r for c, r in recordscalled]

            avgloc1 = loc1.withPos(medianPos1)
            avgloc2 = loc2.withPos(medianPos2)
            ref, alt = bkptRefAltFromPair(avgloc1, avgloc2)
            vcfrec = outvcf.new_record()
            vcfrec.contig = avgloc1.__chrom__
            vcfrec.start = avgloc1.__pos__
            vcfrec.ref = ref
            vcfrec.alts = [alt]
            vcfrec.filter.add(filterstring)
            for key, val in sorted(infoString(callers, make_info_dict(records, medianPos1, medianPos2)).items()):
                vcfrec.info.__setitem__(key, val)
            print(vcfrec, end="", file=outfile)

            if debug:
                for caller, rec in recordscalled:
                    print("#" + str(rec).rstrip() + " ("+caller+")", file=outfile)

    outvcf.close()

def readMergedCalls(infile, filterByChromosome=True, readINFO=False, skipcallers=None):
    """Read a merged callset, and return:
        - dictionary: caller name -> caller idx
        - callsets(list of lists): [calleridx][callidx]
        - calls: callidx -> record from merged"""
    invcf = pysam.VariantFile(infile)
    callerIdx = 0
    callIdx = 0
    callsets = []
    callIdxToCall = []
    callerIdxDict = {}

    if skipcallers is None:
        skipcallers = []

    for rec in invcf.fetch():
        ncalledthis = 0
        if filterByChromosome and not mapped_to_chromosome(rec.chrom):
            continue
        callers = [c for c in rec.info['Callers'] if not c in skipcallers]
        called = []
        for caller in callers:
            if not (caller in called) and not (caller in skipcallers):
                called.append(caller)

                if not caller in callerIdxDict:
                    callerIdxDict[caller] = callerIdx
                    callerIdx += 1
                    callsets.append([])

                callsets[callerIdxDict[caller]].append(callIdx)
                ncalledthis += 1

        assert len(called) == ncalledthis
        if ncalledthis > 0:
            chrom = rec.chrom
            posstart = rec.pos
            callIdxToCall.append((len(called), chrom, posstart, str(rec.REF), str(list(rec.alts)[0]), ",".join(called)))
            callIdx += 1

    return callerIdxDict, callsets, callIdxToCall