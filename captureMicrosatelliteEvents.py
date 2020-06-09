# =============================================================================
# Project : Pindel 
# Py Name: 
# Author :
# Date : 20-05-06
# Email : pengjia@stu.xjtu.edu.cn
# Description : 'Capture complex indels at microsatellite region based on pindel output'
# =============================================================================
import os
import re
import argparse

global file_D, fileDI, DI_ID, D_ID, DIRecords


def getIndelLocation(ref):
    """
    Get the location of indel, if the indel is located in microsatellite region, capture the complex events.
    :param ref: reference  sequence of pindel record.
    :return:
    """
    p = re.compile("[a-z]+")
    delrange = [0, 0]
    for delitem in p.finditer(ref):
        delrange[0] = delitem.start()

        delrange[1] = delitem.end()
    # print(delrange)
    # print(ref[0:0])
    delStr = ref[delrange[0]:delrange[1]].upper()
    elementNum = 1
    delLen = len(delStr)
    refLen = len(ref)
    # print("delLen",delLen)
    if delLen < 1:
        return False, "", "", "", ""
    pos = delrange[0] + delLen
    while pos + delLen < refLen:

        if ref[pos:pos + delLen] == delStr:
            pos = pos + delLen
            elementNum += 1
        else:
            break
    rightPos = pos
    pos = delrange[0]
    while pos - delLen > 0:

        if ref[pos - delLen:pos] == delStr:
            pos = pos - delLen
            elementNum += 1
        else:
            break
    leftPos = pos
    if elementNum > 1:
        return True, leftPos, rightPos, delStr, elementNum
    else:
        return False, "", "", "", ""


def findMismatch(ref, reads, leftPos, rightPos, minAllefrc=0.2):
    """
    Detect if there are substitution of
    :param ref: reference of pindel output.
    :param reads: reads of pindel output.
    :param leftPos: microsatellite left position.
    :param rightPos:microsatellite left position.
    :param minAllefrc: minimum allele fraction.
    :return:
    """
    pos = leftPos - 1
    leftStr = ""
    while pos >= 0:
        # print(pos)
        posBases = [str(read[pos:leftPos]) for read in reads]
        posBases_dict = {}
        allSuportNum = 0
        for base in posBases:
            if base.replace(" ", "") == "" or len(base) < len(leftStr) + 1:
                continue
            allSuportNum += 1
            if base in posBases_dict:
                posBases_dict[base] = posBases_dict[base] + 1
            else:
                posBases_dict[base] = 1
        Mismatch = False
        for base in posBases_dict:
            if posBases_dict[base] / allSuportNum > minAllefrc and base[0] != ref[pos]:
                leftStr = base[0] + leftStr
                Mismatch = True
                continue
        if not Mismatch:
            break
        pos = pos - 1
    pos = rightPos
    rightStr = ""
    while pos < len(ref) - 1:
        pos += 1
        posBases = [str(read[rightPos:pos]) for read in reads]
        posBases_dict = {}
        allSuportNum = 0
        for base in posBases:
            if base.replace(" ", "") == "" or len(base) < len(rightStr) + 1: continue
            allSuportNum += 1
            if base in posBases_dict:
                posBases_dict[base] = posBases_dict[base] + 1
            else:
                posBases_dict[base] = 1

        Mismatch = False
        for base in posBases_dict:
            if posBases_dict[base] / allSuportNum > minAllefrc and base[-1] != ref[pos - 1]:
                rightStr = rightStr + base[-1]
                Mismatch = True
                continue
        if not Mismatch:
            break
    # print(leftStr,rightStr)
    if len(leftStr) < 1 and len(rightStr) < 1:
        return False, "", ""
    else:
        if len(leftStr) >= len(rightStr):
            return True, "left", leftStr
        else:
            return True, "right", rightStr


def processingRecord(RecordInfo, recordNum):
    """
    Processing one deletion record in pindel output
    :param RecordInfo: record string
    :param recordNum: record ID
    :return:
    """
    global file_D, fileDI, DI_ID, D_ID
    print("[Info] Processing Record ", recordNum)
    header = RecordInfo[1]
    headerinfo = header.split("\t")
    ref = RecordInfo[2]
    indelLocation = getIndelLocation(ref)

    if indelLocation[0]:
        BPStart = int(header.split("\t")[4].split(" ")[1])
        # BPEnd = int(header.split("\t")[5])
        BPRagneStart = int(header.split("\t")[6].split(" ")[1])
        # print("HHHH",header.split("\t")[7])
        BPRagneEnd = int(header.split("\t")[7])
        _, leftPos, rightPos, element, elementNum = indelLocation
        reads = [readinfo.split("\t")[0] for readinfo in RecordInfo[3:]]
        label, direction, NT = findMismatch(ref, reads, leftPos=leftPos, rightPos=rightPos)
        # print(elementNum,label, direction, NT)
        if label:
            print("[Info] This is a deletion-insertion event!")
            ref = ref.upper()
            if direction == "left":
                BPStart = BPStart - len(NT)
                BPEnd = BPStart + len(element) + 2
                BPRagneStart = BPRagneStart - len(NT)
                BPRagneEnd = BPRagneEnd

            else:
                RecordInfo[2] = ref[0:leftPos] \
                                + element * (elementNum - 1) + element.lower() \
                                + ref[leftPos + elementNum * len(element) + elementNum:]
                reads = RecordInfo[3:]
                newReads = []
                for read in reads:
                    newReads.append(
                        read[0:leftPos] \
                        + element * (elementNum - 1) + " " * len(element) \
                        + read[leftPos + elementNum * len(element) + elementNum:]
                    )
                # print(RecordInfo[2])
                BPStart = BPStart + elementNum * len(element) - 1
                BPEnd = BPStart + len(NT) + 2
                BPRagneStart = BPRagneStart
                BPRagneEnd = BPEnd
            SVsize = BPEnd - BPStart - 1

            headerinfo[0] = str(DI_ID)
            DI_ID = DI_ID + 1
            headerinfo[1] = "DI " + str(SVsize)
            headerinfo[2] = "NT " + str(len(NT)) + " " + NT
            headerinfo[4] = "BP " + str(BPStart)
            headerinfo[5] = str(BPEnd)
            headerinfo[6] = "BP_range " + str(BPRagneStart)
            headerinfo[7] = str(BPRagneEnd)
            ChrID = headerinfo[3].split(" ")[1]

            if ChrID not in DIRecords:
                DIRecords[ChrID] = {}
            DIRecords[ChrID][BPStart] = RecordInfo[0] + \
                                        "\t".join(headerinfo) + \
                                        RecordInfo[2] + \
                                        "".join(RecordInfo[3])
        else:
            headerinfo[0] = str(D_ID)
            D_ID += 1
            file_D.write(RecordInfo[0])
            file_D.write("\t".join(headerinfo))
            file_D.write(RecordInfo[2])
            file_D.write("".join(RecordInfo[3:]))
    else:
        headerinfo = header.split("\t")
        headerinfo[0] = str(D_ID)
        D_ID += 1
        file_D.write(RecordInfo[0])
        file_D.write("\t".join(headerinfo))
        file_D.write(RecordInfo[2])
        file_D.write("".join(RecordInfo[3:]))


if __name__ == "__main__":

    global file_D, fileDI, D_ID, DI_ID, DIRecords

    ############# read the arguments #################
    parser = argparse.ArgumentParser(
        description='Capture complex indels at microsatellite region based on pindel output')
    parser.add_argument('--pindel_prefix', required=True, type=str, nargs=1,
                        help="path of reference file [required]")
    args = parser.parse_args()
    path_pindel_output_prefix = args.pindel_prefix[0]
    path_output_D = path_pindel_output_prefix + "_D_Simple"
    path_output_DI = path_pindel_output_prefix + "_D_Deletion_Insertion"
    path_output_D_tmp = path_pindel_output_prefix + "_D"
    # path_output_DI_tmp = path_pindel_output_prefix + "_DI"

    #############  initialization #################
    file_D = open(path_output_D, "w")
    file_DI = open(path_output_DI, "w")
    tmpRecordInfo = []
    recordNum = 0
    D_ID = 0
    DI_ID = 0
    DIRecords = {}

    ############# Processing every record in pindel deletion  #################
    for line in open(path_output_D_tmp):
        if "####" in line:

            if len(tmpRecordInfo) == 0:
                tmpRecordInfo.append(line)
            else:
                recordNum += 1
                processingRecord(tmpRecordInfo, recordNum)
                tmpRecordInfo = []
                tmpRecordInfo.append(line)
        elif "ChrID" in line:
            tmpRecordInfo.append(line)
        else:
            tmpRecordInfo.append(line)
    recordNum += 1
    processingRecord(tmpRecordInfo, recordNum)
    file_D.close()

    ############# sort and write to the file  #################

    chrIDList = sorted(list(DIRecords.keys()))
    num = 0
    for chrId in chrIDList:
        thisDIPosList = sorted(list(DIRecords[chrId].keys()))
        for thisDIPos in thisDIPosList:
            lines = DIRecords[chrId][thisDIPos].split("\n")
            headerinfo = lines[1].split("\t")
            headerinfo[0] = str(num)
            num += 1
            lines[1] = "\t".join(headerinfo)
            file_DI.write("\n".join(lines))
    print("[Info] Total capture", DI_ID, "complex events in deletion call set!")
    file_DI.close()
