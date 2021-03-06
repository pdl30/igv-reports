import os
import sys
import json
import math
import argparse
from urllib.request import urlopen
from igv_reports import fasta, ideogram, datauri, tracks, feature, bam, vcf, utils
from igv_reports.varianttable import VariantTable
from igv_reports.bedtable import BedTable


def create_report(args):

    variants_file = args.sites

    if variants_file.endswith(".vcf") or variants_file.endswith (".vcf.gz"):
        table = VariantTable(variants_file, args.goshg2p_url, args.goshg2p_port, args.panel_file,
                             args.info_columns, args.info_columns_prefixes, args.sample_columns,
                             args.lynch_sample)

    elif variants_file.endswith(".bed") or variants_file.endswith(".bed.gz"):
        table = BedTable(variants_file)

    table_json = table.to_JSON()

    session_dict = {}

    # Create file readers for tracks.  This is done outside the loop so initialization happens onc
    trackreaders = []
    if args.tracks is not None:
        for track in args.tracks:
            reader = utils.getreader(track)
            trackreaders.append({
                "track": track,
                "reader": reader
            })


    # loop through variants creating an igv.js session for each one
    for tuple in table.features:

        feature = tuple[0]
        unique_id = tuple[1]

        # Define a genomic region around the variant
        chr = feature.chr
        position = int(math.floor((feature.start + feature.end) / 2)) + 1   # center of region in 1-based coordinates
        start = int (math.floor(feature.start - int(args.flanking) / 2))
        end = int (math.ceil(feature.end + int(args.flanking) / 2))
        region = {
            "chr": chr,
            "start": start,
            "end": end
        }

        # Fasta
        data = fasta.get_data(args.fasta, region)
        fa = '>' + chr + ':' + str(start) + '-' + str(end) + '\n' + data
        fasta_uri = datauri.get_data_uri(fa)
        fastaJson = {
            "fastaURL": fasta_uri,
        }

        # Ideogram
        if(args.ideogram):
            ideo_string = ideogram.get_data(args.ideogram, region)
            ideo_uri = datauri.get_data_uri(ideo_string)
            fastaJson["cytobandURL"] = ideo_uri


        # Initial locus, +/- 20 bases
        initial_locus = chr + ":" + str(position - 20) + "-" + str(position + 20)
        session_json = {
            "locus": initial_locus,
            "reference": fastaJson,
            "tracks": []
        }

        for tr in trackreaders:

            track = tr["track"]
            reader = tr["reader"]
            trackObj = tracks.get_track_json_dict(track)
            data = reader.slice(region)
            trackObj["url"] = datauri.get_data_uri(data)
            if(trackObj["type"] == "alignment"):
                trackObj["height"] = 500

                # Sort TODO -- do this only for SNV
                # if (trackObj["type"]) == "alignment":
                #     trackObj["sort"] = {
                #         "option": "NUCLEOTIDE",
                #         "locus": chr + ":" + str(variant.pos - 1)
                #     }

            session_json["tracks"].append(trackObj)


        # Build the session data URI

        session_string = json.dumps(session_json);

        session_uri = datauri.get_data_uri(session_string)

        session_dict[str(unique_id)] = session_uri

    session_dict = json.dumps(session_dict)

    template_file = args.template
    if None == template_file:
        template_file = os.path.dirname(sys.modules['igv_reports'].__file__) + '/templates/variant_template.html'

    output_file = args.output

    standalone = args.standalone
    with open(template_file, "r") as f:
        data = f.readlines()

        with open(output_file, "w") as o:

            for i, line in enumerate(data):

                if standalone:
                    if line.strip().startswith("<script") and ".js\"" in line:
                        inline_script(line, o, "js", args)
                        continue
                    elif line.strip().startswith("<link") and line.strip().endswith("css\">"):
                        inline_script(line, o, "css", args)
                        continue
                j = line.find('"@TABLE_JSON@"')
                if j >= 0:
                    line = line.replace('"@TABLE_JSON@"', table_json)

                j = line.find('"@SESSION_DICTIONARY@"')
                if j >= 0:
                    line = line.replace('"@SESSION_DICTIONARY@"', session_dict)

                o.write(line)


def inline_script(line, o, source_type, args):
    #<script type="text/javascript" src="https://igv.org/web/test/dist/igv.min.js"></script>
    if source_type == "js":
        s = line.find('src="')
        offset = 5
        o.write('<script type="text/javascript">\n')
    elif source_type == "css":
        s = line.find('href="')
        offset = 6
        o.write('<style type="text/css">\n')
    else:
        raise KeyError("Inline script must be either js- or css-file")
    if s > 0:
        e = line.find('">', s)
        # url = line[s+offset:e]
        js_file = open(args.igv_js, 'r')
        content = js_file.read()
        js_file.close()
        o.write(content)
        if source_type == "js":
            o.write('</script>\n')
        else:
            o.write('</style>\n')
    else:
        raise ValueError("No file path in {l} for inline script.".format(l=line))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sites", help="vcf file defining variants, required")
    parser.add_argument("fasta", help="reference fasta file, required")
    parser.add_argument("--ideogram", help="ideogram file in UCSC cytoIdeo format")
    parser.add_argument("--tracks", nargs="+", help="list of track files")
    parser.add_argument("--template", help="html template file", default=None)
    parser.add_argument("--output", help="output file name", default="igvjs_viewer.html")
    parser.add_argument("--info-columns", nargs="+", help="list of VCF info field names to include in variant table")
    parser.add_argument("--info-columns-prefixes", nargs="+", help="list of prefixes of VCF info field names to include in variant table")
    parser.add_argument("--sample-columns", nargs="+", help="list of VCF sample/format field names to include in variant table")
    parser.add_argument("--flanking", help="genomic region to include either side of variant", default=1000)
    parser.add_argument('--standalone', help='Print more data', action='store_true')
    parser.add_argument('--igv_js', help='Path of igv.min.js file')
    parser.add_argument('--lynch_sample', help='Does the sample have the LYNCH panel applied', action='store_true')
    parser.add_argument('--panel_file', help='Bed file for panel; file should have gene name on column 4')
    parser.add_argument('--goshg2p_url', help='GOSHG2P URL e.g. 10.101.XX.XX', )
    parser.add_argument('--goshg2p_port', help='GOSHG2P PORT e.g. 3306')
    args = parser.parse_args()
    create_report(args)

if __name__ == "__main__":

    main()
