import io
import json
import pysam
from .feature import Feature
import requests


class VariantTable:

    # Always remember the *self* argument
    def __init__(self, vcfFile, panel_file=True, info_columns = None, info_columns_prefixes = None, sample_columns = None,
                 lynch_sample=False):

        vcf = pysam.VariantFile(vcfFile)
        self.panel_file = panel_file
        self.info_fields =  info_columns or []
        self.info_field_prefixes = info_columns_prefixes or []
        self.sample_fields = sample_columns or []
        self.variants = []
        self.features = []   #Bed-like features
        self.genes = self.load_genes()
        self.lynch_sample = lynch_sample
        for record in vcf.header.records:
            if 'ALAMUT_ANN' in str(record):
                alamut_ann = str(record).rstrip()
                self.alamut_ann = alamut_ann.split('Format: ')[1].split('|')

        lynch_variant_added = False
        variant_id = 0
        for unique_id, var in enumerate(vcf.fetch()):
            try:
                if lynch_sample and not lynch_variant_added and int(var.chrom.strip('chr')) > 1:
                    if (var.chrom == 'chr2' and var.pos > 47641560) or int(var.chrom.strip('chr')) >= 3:
                        # You should now add the lynch_variant
                        self.features.append((Feature('chr2', 47641559, 47641561, ''), variant_id))
                        variant_id += 1
                        lynch_variant_added = True
            except ValueError as e:
                pass

            self.variants.append((var, variant_id))
            chr = var.chrom
            start = var.pos - 1
            end = start + 1       #TODO -- handle structure variants and deletions > 1 base
            self.features.append((Feature(chr, start, end, ''), variant_id))
            variant_id += 1

    def to_JSON(self):
        # In order to create correct HGVS, I need to know the correct gene,
        # For correct gene, you should have this in the panel file
        # Longest transcript should be in VCF
        # Selected transcript will need an API
        json_array = [];
        lynch_variant_added = False
        variant_id = 0
        for variant, unique_id in self.variants:
            gene_name = self.get_gene_name(variant)
            gene_transcripts = self.goshg2p_connection(gene_name)

            if gene_transcripts:
                selected_transcript = gene_transcripts[0]['transcript']
            else:
                selected_transcript = self.get_longest_transcript(variant)
            hgvs = self.get_hgvs(variant, selected_transcript)
            obj = {
                'unique_id': variant_id,
                'CHROM': variant.chrom,
                'POSITION': variant.pos,
                'REF': variant.ref,
                'ALT': ','.join(variant.alts),
                'ID': hgvs,
                'GENE': gene_name
            }

            if variant.id is not None:
                obj['ID'] = render_ids(variant.id)

            for h in self.info_fields:
                v = ''
                if h in variant.info:
                    if h == 'ANN':
                        genes, effects, impacts, transcript, gene_id, aa_alt, nt_alt = decode_ann(variant)
                    elif h == 'COSMIC_ID':
                        v = render_id(v)
                    else:
                        v = render_values(variant.info[h])
                if h == 'ANN':
                    obj['GENE'] = genes
                    obj['EFFECTS'] = effects
                    obj['IMPACT'] = impacts
                    obj['TRANSCRIPT'] = transcript
                    obj['GENE_ID'] = gene_id
                    obj['PROTEIN ALTERATION'] = aa_alt
                    obj['DNA ALTERATION'] = nt_alt
                else:
                    obj[h] = v

            for h in self.info_field_prefixes:
                v = ''
                for field in variant.info:
                    if field.startswith(h):
                        obj[field] = render_values(variant.info[field])

            for h in self.sample_fields:
                for sample, values in variant.samples.items():
                    v = ''
                    try:
                        v = values[h]
                    except KeyError:
                        # ignore if key is not present
                        pass

                    obj[f'{sample}:{h}'] = render_values(v)
            try:
                if self.lynch_sample and not lynch_variant_added and int(variant.chrom.strip('chr')) > 1:
                    if (variant.chrom == 'chr2' and variant.pos > 47641560) or int(variant.chrom.strip('chr')) >= 3:
                        # You should now add the lynch_variant
                        json_array.append({'unique_id': variant_id, 'CHROM': 'chr2', 'POSITION': 47641560,
                                           'REF': 'A', 'ALT': 'T', 'ID': 'LYNCH FOR REVIEW - variant not called',
                                           'GENE': 'MSH2'})
                        variant_id += 1
                        obj['unique_id'] = variant_id
                        lynch_variant_added = True
            except ValueError as e:
                pass
            json_array.append(obj)
            variant_id += 1

        if not any(obj['ID'] for obj in json_array):
            # Remove ID column if none of the records actually had an ID.
            for obj in json_array:
                del obj['ID']

        return json.dumps(json_array)

    def load_genes(self):
        genes = []
        with open(self.panel_file) as f:
            for line in f:
                word = line.rstrip().split('\t')
                if word[3] not in genes:
                    genes.append(word[3])
        return genes

    def get_gene_name(self, variant):
        for ann in variant.info['ALAMUT_ANN']:
            zipped = zip(self.alamut_ann, ann.split('|'))
            format_dict = {f[0]: f[1] for f in zipped}
            if 'gene' in format_dict:
                if format_dict['gene'] in self.genes:
                    return format_dict['gene']
        return None  # If no gene, don't bother trying to annotate HGVS

    def get_longest_transcript(self, variant):
        """
        Loops over the annotations and a dict with the longest transcript info
        :param variant:
        :return:
        """
        longest_transcript = {'transcript': '', 'length': 0}
        for ann in variant.info['ALAMUT_ANN']:
            zipped = zip(self.alamut_ann, ann.split('|'))
            format_dict = {f[0]: f[1] for f in zipped}
            if 'transLen' in format_dict:
                if int(format_dict['transLen']) > int(longest_transcript['length']):
                    longest_transcript['transcript'] = format_dict['transcript']
                    longest_transcript['length'] = int(format_dict['transLen'])
        if longest_transcript:
            return longest_transcript['transcript']
        else:
            return

    def get_hgvs(self, variant, transcript):
        for ann in variant.info['ALAMUT_ANN']:
            zipped = zip(self.alamut_ann, ann.split('|'))
            format_dict = {f[0]: f[1] for f in zipped}
            if 'transcript' in format_dict:
                if transcript in format_dict['transcript']:
                    return format_dict['cNomen']

    @staticmethod
    def goshg2p_connection(gene_name):
        page = 1
        last_page = False
        all_results = []
        while not last_page:
            json_poll_success = False
            while not json_poll_success:
                MAX_RETRIES = 20
                session = requests.Session()
                adapter = requests.adapters.HTTPAdapter(max_retries=MAX_RETRIES)
                session.mount("http://", adapter)
                response = session.get(
                    url=f"http://10.101.45.28:8014/goshg2p/api/transcripts/{gene_name}?page={page}")
                try:
                    response_json = response.json()
                    response_status = response.status_code
                    if response_status != 200:
                        raise ValueError(f'status code = {response_status}')
                    json_poll_success = True
                    all_results.extend(response_json['results'])
                except (json.JSONDecodeError, ValueError) as e:
                    print(e)
                    raise
            if response_json["next"]:
                page += 1
            else:
                last_page = True
        return all_results

def render_value(v):
    """Render given value to string."""
    if isinstance(v, float):
        # ensure that we don't waste space by insignificant digits
        return f'{v:.2g}'
    elif v.startswith('http://') or v.startswith('https://'):
        return create_link(v)
    else:
        return str(v)


def render_values(v):
    if isinstance(v, str) or isinstance(v, int) or isinstance(v, float):
        return render_value(v)
    return ','.join(map(render_value, v))


def render_id(v):
    if v.startswith('COSM'):
        return (f'<a href = "https://cancer.sanger.ac.uk/cosmic/mutation/overview?'
                f'id={v[4:]}" target="_blank">{v}</a>')
    return v


def render_ids(v):
    return ','.join(map(render_id, v.split(';')))


def decode_ann(variant):
    """Decode the standardized ANN field to something human readable."""
    annotations = ([variant.info['ANN'].split('|'
                   )] if isinstance(variant.info['ANN'],
                   str) else [e.split('|') for e in variant.info['ANN']])
    genes = []
    effects = []
    impacts = []
    transcripts = []
    gene_ids = []
    aa_alts = []
    nt_alts = []
    for allele in variant.alts:
        for ann in annotations:
            ann_allele, kind, impact, gene, gene_id = ann[:5]
            feature_id = ann[6]
            nt_mod, aa_mod = ann[9:11]

            if allele != ann_allele:
                continue

            full = '|'.join(ann)
            # Keep the most severe effect.
            # Link out to Genecards and show the full record in a tooltip.
            genes.append(gene)
            gene_ids.append(gene_id)
            effects.append(kind.replace('&', '/'))
            impacts.append(impact)
            transcripts.append(feature_id)
            aa_alts.append(aa_mod)
            nt_alts.append(nt_mod)
    return ','.join(genes), ','.join(effects), ','.join(impacts), ','.join(transcripts), ','.join(gene_ids), ','.join(aa_alts), ','.join(nt_alts)

def create_link(url):
    """Create an html link for the given url"""
    return (f'<a href = "{url}" target="_blank">{url}</a>')