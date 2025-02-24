#!/usr/bin/env python3
# Author:  Rachel Erhlich, edited by Sam Czerski for Microbiome data

from glob import iglob, glob
from subprocess import Popen, PIPE, call, check_output, DEVNULL
import os
import sys
import re
from collections import defaultdict, Counter, namedtuple, OrderedDict
import argparse
from time import sleep
import json

import pandas as pd
import numpy as np


import plotly.io
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.backends.backend_pdf as plt_pdf
from csv2pdf import convert


#sam required packages
sys.path.append('/home/sam/miniconda3/pkgs')
sys.path.append('/home/bin/python_modules')
sys.path.append('/home/sam/miniconda3/lib/python3.9/site-packages')


class JsonToReports:
    def __init__(self, report_pbcromwell_dir, outdir):
        self.report_pbcromwell_dir = report_pbcromwell_dir
        self.outdir = outdir
        self.data = None

        self.json_path = self._get_json_path()
        self.nickname = self.json_path.split('/')[-1].split('.')[-3]

        self._get_data()
        self.tables = self._get_tables()
        self.images = self._get_images()




    def _get_json_path(self):
        json_path = glob(self.report_pbcromwell_dir + '/*.report.json')
        assert len(json_path) == 1
        json_path = json_path[0]
        return json_path

    def _get_data(self):
        with open(self.json_path) as f:
            self.data = json.load(f)

    def _get_tables(self):
        metrics = dict()
        for attribute in self.data["attributes"]:
            metrics[attribute['name']] = attribute['value']
        #metrics['data'] = self.report_pbcromwell_dir
        metrics = pd.DataFrame(pd.Series(metrics))#.T
        metrics.index.name = 'Metric'
        metrics.reset_index(inplace=True)
        metrics.columns = ['Metric', 'Value']

        tables = dict()
        tables[self.nickname + ' summary'] = metrics

        for json_table in self.data['tables']:
            table = dict()
            for column in json_table['columns']:
                table[column['header']] = column['values']
            table = pd.DataFrame(table)


            if table.columns[:2].to_list() == ['Sample Name', 'Barcode']:
                if (table['Sample Name'] == table['Barcode']).head(len(table) - 1).all():
                    assert table.iloc[-1]['Sample Name'] == 'No Name'
                    assert table.iloc[-1]['Barcode'] == 'Not Barcoded'
                    table = table.drop(columns=['Sample Name'])


        return tables

    def _get_images(self):
        non_thumbnail_plots = dict()
        for group_of_plots in self.data['plotGroups']:
            for single_plot in group_of_plots['plots']:

                if single_plot['plotType'] == 'image':
                    non_thumbnail_plots[single_plot['caption']] = self.report_pbcromwell_dir + '/' + single_plot['image']
                elif single_plot['plotType'] == 'plotly':
                    image_json_path = single_plot['image']
                    assert image_json_path.endswith('.json.gz')
                    image_json_path = self.report_pbcromwell_dir + '/' + image_json_path.removesuffix('.gz')
                    if not os.path.isfile(image_json_path):
                        cmd = ['gunzip', '-f', image_json_path + '.gz']
                        call(cmd)
                    assert os.path.isfile(image_json_path)

                    image_png_path = image_json_path.removesuffix('.json') + '.png'
                    if not os.path.isfile(image_png_path):
                        fig = plotly.io.read_json(image_json_path)
                        plotly.io.write_image(fig, image_png_path)
                    assert os.path.isfile(image_png_path)
                    non_thumbnail_plots[single_plot['caption']] = image_png_path
        return non_thumbnail_plots

    def save_to_pdf(self):
        with plt_pdf.PdfPages(self.outdir + '/' + self.nickname + '_plots.pdf') as pdf:
            for table_name, table in self.tables.items():
                fig = plt.figure()
                plt.axis('off')
                tab = plt.table(cellText=table.values, colLabels=table.columns, loc='center')
                tab.auto_set_font_size(False)
                tab.set_fontsize(3) # 3 is good for the demux tables
                plt.title(table_name)
                pdf.savefig(fig)
                plt.close()

            for caption, png_path in self.images.items():
                img = mpimg.imread(png_path)
                fig = plt.figure()
                plt.imshow(img)
              #  plt.title(caption) unnecessary for these plots, but I think I need it for others
                plt.axis('off')
                pdf.savefig(fig)
                plt.close()

    def save_to_excel(self):
        print(self.outdir)
        print(self.nickname)
        with pd.ExcelWriter(self.outdir + '/' + self.nickname + '_results.xlsx', engine='xlsxwriter') as writer:
            for table_name, table in self.tables.items():
                table.to_csv(self.outdir + '/' + table_name + '.csv', index=False)
                table.to_excel(writer, sheet_name=table_name, index=False)

                sheet = writer.sheets[table_name]
                sheet.set_zoom(200)
                sheet.autofit()

            worksheet = writer.book.add_worksheet('plots')
            for i, (caption, png_path) in enumerate(self.images.items(), 1):
                worksheet.insert_image('A' + str((i * 40) - 39), png_path)


    def get_metrics(self, cell_path, job_num):
        barcode_summary_path = glob(self.outdir + '/outputs/*barcodes_summary.csv')
        assert len(barcode_summary_path) == 1
        barcode_summary_path = barcode_summary_path[0]

        metrics = pd.read_csv(barcode_summary_path)

        metrics.drop(columns=['Sample Name'], inplace=True)
        metrics['library_strategy'] = 'WGS'
        metrics['library_source'] = 'GENOMIC'
        metrics['library_selection'] = 'RANDOM'
        metrics['library_layout'] = 'single'
        metrics['platform'] = 'PACBIO_SMRT'

        metrics['filetype'] = 'bam'
        metrics['assembly'] = 'unaligned' # submitting a bam w/o a ref based assembly

        metrics['cell_path'] = cell_path
        metrics['demux_dir'] = self.outdir

        hifi_xml = metrics.apply(lambda x: self.get_single_barcode_xml_path(x), axis=1)
        metrics['hifi_xml'] = hifi_xml

        metrics['hifi_bam'] = metrics['hifi_xml'].str.replace('consensusreadset.xml', 'bam')
        metrics['filename'] = metrics['hifi_bam']
        for i, row in metrics.iterrows():
            if row['Barcode'] != 'Not Barcoded':
                assert os.path.isfile(row['hifi_bam'])

        metrics['barcode_xml'] = self.args.barcode_xml

        metrics['demux_software_version'] = self.get_demux_software_version(metrics['hifi_xml'].iloc[0])

        assert '/data/pacbio/sequel_IIe/' in metrics['cell_path'].iloc[0]
        metrics['instrument_model'] = 'PacBio Sequel IIe'

        metrics.to_csv(self.outdir + '/results_summary.csv', index=False)
        metrics.to_csv('demux_result_' + job_num + '.txt', index=False)

    def get_demux_software_version(self, hifi_xml_path):
        with open(hifi_xml_path, 'r') as f:
            for line in f:
                if '<pbmeta:VersionInfo Name="smrtlink" Version=' in line:
                    version = line.split('"')[-2]
                    return version
        raise NotImplementedError

    @staticmethod
    def get_single_barcode_xml_path(row):
        xml = glob(row['demux_dir'] + '/outputs/*.' + row['Barcode'] + '.consensusreadset.xml')

        if len(xml) == 0:
            assert row['Barcode'] == 'Not Barcoded'
            return ''

        assert len(xml) == 1
        xml = xml[0]

        xml = os.path.realpath(xml)
        assert os.path.isfile(xml)

        return xml

    @staticmethod
    def per_bc_report_pdf(barcodes_summary_csv, outdir):
        #FIX LATER
        with plt_pdf.PdfPages(outdir + '/' + 'barcodes_summary.pdf') as pdf:
            fig = plt.figure(figsize=(8.27,11))
            data = pd.read_csv(barcodes_summary_csv)
            data = pd.DataFrame(data)
            tab = plt.table(cellText=data.values, colLabels=data.columns, loc='center', colWidths=[0.1 for col in range(data.columns.size)])
            plt.axis('off')
            tab.auto_set_font_size(False)
            tab.set_fontsize(3)
            plt.title("barcodes summary")
            pdf.savefig(fig)
            plt.close()
        
    @staticmethod
    def convert_csv2pdf(barcodes_summary_csv, outdir):
        
        pre_csv = pd.read_csv(barcodes_summary_csv)
        post_csv = pre_csv.apply(lambda x: pd.Series(x.dropna().values))

        cols_float_to_int = ['Barcode Quality','HiFi Reads','HiFi Read Length (mean, bp)','HiFi Yield (bp)', 'Polymerase Read Length (mean, bp)', 'Polymerase Yield (bp)']
        for c in cols_float_to_int:
            post_csv[c] = post_csv[c].apply(np.int64)

        post_csv.to_csv(outdir + '/' + 'per_bc_report.csv', sep=',', index=False)

        convert(outdir + '/' + "per_bc_report.csv", outdir + '/' + 'barcodes_summary.pdf', orientation="L")

def make_args():
    parser = argparse.ArgumentParser(description='',
                                     add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    required_flags = parser.add_argument_group('Required arguments')

    required_flags.add_argument('-c', '--ccs_and_demux_dir',
                                help='',
                                required=True)

    required_flags.add_argument('-s', '--sequel_platform',
                                help='',
                                required=True)

    required_flags.add_argument('-o', '--report_output_dir',
                                help='',
                                required=True)

    return parser.parse_args()


def main():
    args = make_args()
    sequel_platform="Sequel_I"

    ccs_demux_outdir = args.ccs_and_demux_dir
    report_outdir = args.report_output_dir

    data_dir = glob(ccs_demux_outdir + '/ccs' + '/cromwell-executions/pb_ccs/*/call-pbreports_ccs2/execution')

    data_dir = data_dir[0]
        
    ccs_report = JsonToReports(data_dir, report_outdir)
    ccs_report.save_to_excel()
    ccs_report.save_to_pdf()

    if sequel_platform == "Sequel_I":
        report_dir = glob(ccs_demux_outdir + '/demux_no_peek' + '/cromwell-executions/pb_demux_ccs/*/call-demultiplex_barcodes/demultiplex_barcodes/*/call-barcode_report/execution')
        assert len(report_dir) == 1
        report_dir = report_dir[0]
    else:
        report_dir = glob(ccs_demux_outdir + '/demux_no_peek' + '/cromwell-executions/pb_demux_ccs/*/call-barcode_report/execution')
        assert len(report_dir) == 1
        report_dir = report_dir[0]

    demux_report = JsonToReports(report_dir, report_outdir)
    demux_report.save_to_excel()
    demux_report.save_to_pdf()

    if sequel_platform == "Sequel_I":
        per_bc_summary = glob(ccs_demux_outdir + '/demux_no_peek' + '/outputs/' + '/barcode_ccs_summary.csv')
        assert len(per_bc_summary) == 1
        per_bc_summary = per_bc_summary[0]

    else:
        per_bc_summary = glob(ccs_demux_outdir + '/demux_no_peek' + '/cromwell-executions/pb_demux_ccs/*/call-demultiplex_barcodes/demultiplex_barcodes/*/call-barcode_report/execution' + '/barcode_ccs_summary.csv')
        print(per_bc_summary, "2")
        assert len(per_bc_summary) == 1
        per_bc_summary = per_bc_summary[0]

    #JsonToReports.per_bc_report_pdf(per_bc_summary, report_outdir)
    JsonToReports.convert_csv2pdf(per_bc_summary, report_outdir)
            #job_num = args.subread_path_file.split('_')[-1].split('.')[0]
        #demux_job.get_metrics(cell_dir, job_num)


if __name__ == '__main__':
    main()
