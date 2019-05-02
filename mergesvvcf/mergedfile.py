from __future__ import print_function
import pyximport

pyximport.install()
import os
import operator
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
    if type(val) is list or type(val) is tuple:
        val = val[0]
    try:
        i = int(val)
    except ValueError:
        i = val
    return i

    # def make_info_dict(callers, records, pos1, pos2):
    #     """ generate 'median' results for info fields from all records """
    #     callermap = {
    #         "delly": {"fields": ["DV", "RV"], "tumor": 0},
    #         "svaba": {"fields": ["SR", "DR"], "tumor": 1},
    #         "gridss": {"fields": ["RP", "SR"], "tumor": 1},
    #         "brass": {"fields": ["PS", "RC"], "tumor": 1},
    #         "smoove": {"fields": ["SR", "PE"], "tumor": 0},
    #     }
    #     # fields = ['CHR2', 'END', 'SVTYPE', 'SVLEN'] + ["SR", "DR"] + ["DV", "RV"] + ["RP"]
    #     fields = [x["fields"] for x in callermap.values()]
    #     fields = [item for sublist in fields for item in sublist]
    #     info = {}
    #     for field in fields:
    #         answers = []
    #         for caller, record in zip(callers, records):
    #             if caller in callermap.keys() and field in callermap[caller]["fields"]:
    #                 if field in record.format:
    #                     answers.append([caller,int_if_possible(record.samples[callermap[caller]["tumor"]][field])])
    #                 elif field in record.info:
    #                     answers.append([caller,int_if_possible(record.info[field])])
    #                 # elif field in record.format:
    #                 #     answers.append([caller,int_if_possible(record.samples[callermap[caller]["tumor"]][field])])
    #             nanswers = len(answers)
    #         if nanswers > 0:
    #             # sorted_answers = sorted(answers)
    #             # # doesn't quite deal with even #s correctly - can't average strings
    #             # median_pos = int(nanswers / 2)
    #             # median_answer = sorted_answers[median_pos]
    #             # if not median_answer == 0:
    #             #     info[field] = median_answer
    #             for a in answers:
    #                 info[a[0] + "_" + field] = a[1]

    if "SVTYPE" in info and info["SVTYPE"] in ["DUP", "DUP:TANDEM", "DEL", "INV"]:
        if not "SVLEN" in info:
            info["SVLEN"] = pos2 - pos1
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


def getSVTYPE(chr1, chr2, extend1, extend2):
    """Get SVTYPE from extend right"""
    if chr1 != chr2:
        return "TRA"

    eventmap = {True: {True: "INV", False: "DUP"}, False: {True: "DEL", False: "INV"}}

    return eventmap[extend1][extend2]


def merge(
    filenames,
    programs,
    forceSV,
    outfile,
    slop=0,
    verbose=True,
    output_ncallers=False,
    min_num_callers=0,
    filterByChromosome=True,
    noFilter=False,
    debug=False,
):
    """Merge several VCFs from different programs into a new VCF file."""

    outvcf = pysam.VariantFile(os.devnull, "w")
    if min_num_callers > 0:
        outvcf.header.filters.add(
            "LOWSUPPORT", None, None, "Not called by enough callers in ensemble"
        )

    # Returns true if the variant is PASS in the VCF file
    def passed_variant(record):
        """Did this variant pass?"""
        return (
            "PASS" in list(record.filter) or len(list(record.filter)) == 0 or noFilter
        )

    # def infoString(callers, infodict):
    #     """
    #     Generate an INFO string from the INFO dictionary plus
    #     the list of callers.
    #     """
    #     info = {}
    #     # fields = ['SVTYPE', 'SVLEN']
    #     # fields = ["svaba_SR", "svaba_DR"] + ["delly_DV", "delly_RV"] + ["gridss_RP", "gridss_SR"]
    #     # for field in fields:
    #     for field in infodict.keys():
    #         if field in infodict:
    #             res = infodict[field]
    #             if type(res) is list:
    #                 res = res[0]
    #             info[field] = res
    #     if output_ncallers:
    #         info["NumCallers"] = int(len(set(callers)))
    #     info["Callers"] = ",".join(list(set(callers)))
    #     return info

    samplemap = {}  # the samples in each file
    infomap = {}  # the info fields in each file
    formatmap = {}  # the format fields in each file
    outsamples = set()  # standardized samples

    calldict = variantdict.variantmap(awindow=0, svwindow=slop)
    for (infile, program) in zip(filenames, programs):
        count = 0
        try:
            vcf_reader = pysam.VariantFile(infile)

            # Add contigs from all files
            for contig in sorted(vcf_reader.header.contigs):
                try:
                    outvcf.header.contigs.add(contig)
                except ValueError:
                    pass

            outvcf.header.add_meta(
                key="Caller", items=[("ID", program), ("File", infile)]
            )

            # Add filters from all files
            for hfilter in vcf_reader.header.filters.itervalues():
                if hfilter.name not in ["PASS", "LOWSUPPORT"]:
                    outvcf.header.filters.add(
                        "{}_{}".format(program, hfilter.name),
                        None,
                        None,
                        hfilter.description,
                    )
            # Add info fields from all files
            infos = []
            for info in vcf_reader.header.info.itervalues():
                if info.name not in ["SVTYPE", "SVLEN", "END", "CHR2", "NumCallers"]:
                    outvcf.header.info.add(
                        "{}_{}".format(program, info.name),
                        info.number,
                        info.type,
                        info.description,
                    )
                    infos.append(info.name)
            infomap[program] = infos
            # Add format fields from all files
            formats = []
            for hformat in vcf_reader.header.formats.itervalues():
                num = None
                if hformat.name == "GT" or hformat.number == "G":
                    num = "G"
                    continue  # Genotype is special cant do multiple
                outvcf.header.formats.add(
                    "{}_{}".format(program, hformat.name),
                    num or hformat.number,
                    hformat.type,
                    hformat.description,
                )
                formats.append(hformat.name)
            formatmap[program] = formats
            # Collect samples all files
            samples = {}
            for meta in vcf_reader.header.records:
                if meta.key == "SAMPLE" and "SampleName" in meta.keys():
                    samples[meta["ID"]] = meta["SampleName"]
            if not samples:
                for samp in vcf_reader.header.samples:
                    samples[samp] = samp.rsplit("/")[-1].replace(".bam", "")
            samplemap[program] = samples
            outsamples.update(samples.values())

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
                raise (e)
            else:
                pass

    # Add samples
    for samp in sorted(outsamples):
        outvcf.header.samples.add(samp)

    # Write the results in a master vcf file for the sample
    outvcf.header.info.add("Callers", ".", "String", "Callers that made this call")

    if output_ncallers:
        outvcf.header.info.add(
            "NumCallers", 1, "Integer", "Number of callers that made this call"
        )

    outvcf.header.info.add(
        "CHR2",
        "1",
        "String",
        "Chromosome for END coordinate in case of a translocation",
    )
    outvcf.header.info.add(
        "END", "1", "Integer", "End position of the structural variant"
    )
    outvcf.header.info.add("SVTYPE", "1", "String", "Type of structural variant")
    outvcf.header.info.add("SVLEN", "1", "String", "Length of structural variant")
    outvcf.header.info.add(
        "STRANDS",
        "1",
        "String",
        "Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)",
    )

    print(outvcf.header, end="", file=outfile)

    for variant in sorted(calldict, key=operator.itemgetter(0)):
        callers = variant[2]
        num_callers = len(set(callers))
        passes = num_callers >= min_num_callers
        filterstring = "." if passes else "LOWSUPPORT"

        if len(variant) == 3:  # snv/indel
            loc, allele, callers = variant
            if allele is None:
                print("Allele is none: loc, allele, callers = ", loc, allele, callers)
                continue
            chrom, pos, _, _ = loc.asTuple()
            vcfrec = outvcf.new_record()
            vcfrec.contig = chrom
            vcfrec.pos = pos
            vcfrec.ref = allele[0]
            vcfrec.alts = [allele[1]]
            vcfrec.filter.add(filterstring)
            vcfrec.info["Callers"] = ",".join(list(set(callers)))
            if output_ncallers:
                vcfrec.info["NumCallers"] = str(len(set(callers)))
            print(vcfrec, end="", file=outfile)
        else:
            loc1, loc2, callers, medianPos1, medianPos2, recordscalled = variant
            if filterByChromosome and not mapped_to_chromosome(loc2.chrom):
                continue

            # records = [r for c, r in recordscalled]

            avgloc1 = loc1.withPos(medianPos1)
            avgloc2 = loc2.withPos(medianPos2)
            extend1, extend2 = avgloc1.__right__, avgloc2.__right__
            ref, alt = bkptRefAltFromPair(avgloc1, avgloc2)
            vcfrec = outvcf.new_record()
            vcfrec.contig = avgloc1.chrom
            vcfrec.pos = avgloc1.pos
            vcfrec.ref = ref
            vcfrec.alts = [alt]
            vcfrec.filter.add(filterstring)
            # for key, val in sorted(infoString(callers, make_info_dict(callers, records, medianPos1, medianPos2)).items()):
            #     vcfrec.info.__setitem__(key, val)
            vcfrec.info["Callers"] = ",".join(list(set(callers)))
            if output_ncallers:
                vcfrec.info["NumCallers"] = len(set(callers))
            vcfrec.info["CHR2"] = avgloc2.chrom
            vcfrec.stop = avgloc2.pos
            svtype = getSVTYPE(avgloc1.chrom, avgloc2.chrom, extend1, extend2)
            vcfrec.info["SVTYPE"] = svtype
            if avgloc1.chrom == avgloc2.chrom:
                vcfrec.info["SVLEN"] = str(avgloc2.pos - avgloc1.pos)
            strandmap = {False: "+", True: "-"}
            vcfrec.info["STRANDS"] = strandmap[extend1] + strandmap[extend2]
            for caller, rec in recordscalled:
                if noFilter:
                    for fil in rec.filter:
                        if not fil in [".", "PASS"]:
                            vcfrec.filter.add("{}_{}".format(caller, fil))
                for info in infomap[caller]:
                    try:
                        vcfrec.info["{}_{}".format(caller, info)] = rec.info[info]
                    except KeyError:
                        continue
                    except:
                        print(info, rec.info[info])
                        raise
                for samp in samplemap[caller].items():
                    for form in formatmap[caller]:
                        try:
                            vcfrec.samples[samp[1]][
                                "{}_{}".format(caller, form)
                            ] = rec.samples[samp[0]][form]
                        except KeyError:
                            continue
                        except:
                            print(form, rec.samples[samp[0]][form])
                            raise
            print(vcfrec, end="", file=outfile)
            if debug:
                for caller, rec in recordscalled:
                    print("#" + str(rec).rstrip() + " (" + caller + ")", file=outfile)

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
        callers = [c for c in rec.info["Callers"] if not c in skipcallers]
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
            callIdxToCall.append(
                (
                    len(called),
                    chrom,
                    posstart,
                    str(rec.REF),
                    str(list(rec.alts)[0]),
                    ",".join(called),
                )
            )
            callIdx += 1

    return callerIdxDict, callsets, callIdxToCall
