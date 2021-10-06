# -*- coding: utf-8 -*-
import logging
import os
import sys
import webbrowser
from multiprocessing import freeze_support
from time import sleep

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5 import __file__ as qt_file
from PyQt5.QtWidgets import QMainWindow, QFileDialog

import pageant as main

main.plt.switch_backend('Agg')
if main.platform == 'Darwin':
    os.chdir(main.raw_dir)
default_font_size = 10 if main.platform == 'Linux' else 11
y_offset = 0 if main.platform == 'Windows' else 10
QtGui.QGuiApplication.addLibraryPath(os.path.join(os.path.dirname(qt_file), 'Qt', 'plugins', 'platforms'))
QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)

default_config = {value: key[value] for key in main.load_config().values() for value in key}
default_config.update({
    'name': 'HG001',
    'input_file': './personal_genome/HG001.vcf.gz',
    'config_file': './bin/config.ini',
    'quan_data': './algorithm/Quantitative',
    'qual_data': './algorithm/Qualitative',
    'quan_ref': './population_genome',
    'qual_ref': './population_genome',
    # 'qr_snps_txt': './personal_genome/fingerprint_snps.txt',
    'maf_ref': './personal_genome/g1k.vcf.gz',
    'ps_ref': './personal_genome/g1k.vcf.gz',
    'concord_ref': './personal_genome/concordance.vcf.gz'
})


def index_to_sep(ind: int) -> str:
    return {0: '\t', 1: ',', 2: ' '}[ind]


def umap(*args, **kwargs):
    main.get_plink_dir(os.path.join(os.getcwd(), 'bin'))
    main.ps_analyse(*args, **kwargs)


class GUIHandler(logging.Handler):
    def __init__(self, connector: QtWidgets.QProgressDialog, progress_value, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.connector = connector
        self.progress = progress_value

    def emit(self, record):
        msg = self.format(record)
        if 'Progress of the analysis: ' in msg:
            self.progress.value = float(msg.split(': ')[-1].replace('%', ''))
        elif 'start' in msg:
            self.connector.setLabelText(msg)


class MyThread(QtCore.QThread):
    def __init__(self, handler, parameters: dict):
        super().__init__()
        self.statue = False
        self.res = 'Analysis failed to start.'
        self.parameters = parameters
        self.handler = handler
        self.setStackSize(10240000)
        main.gui_log_setup(self.handler)

    def run(self):
        sleep(0.5)
        try:
            self.res = main.main(**self.parameters)
        except Exception as e:
            self.res = str(e).capitalize()
        else:
            self.statue = True
        finally:
            self.quit()
            main.gui_log_remove(self.handler)


class API_Thread(QtCore.QThread):
    def __init__(self, func: main.Callable, *args,
                 end_: str = 'Finish! Result has been saved in the working directory.', **kwargs):
        super().__init__()
        self.statue = False
        self.res = 'Analysis failed to start.'
        self.func = func
        self.args = args
        self.kwargs = kwargs
        self.end = end_
        self.setStackSize(10240000)

    def run(self):
        sleep(0.5)
        try:
            self.func(*self.args, **self.kwargs)
        except Exception as e:
            self.res = str(e).capitalize()
        else:
            self.statue = True
            self.res = self.end
        finally:
            self.quit()


class ProgressValue:
    def __init__(self):
        self.value = 0


class ProgressNow:
    def __init__(self, inter: int):
        self.now_value = 0
        self.inter = 1000 / inter

    def value(self):
        return int(self.now_value // self.inter)

    def add(self):
        self.now_value += 1

    def set(self, value: int):
        self.now_value = value * self.inter


class Ui_PAGEANT(object):
    def setupUi(self, PAGEANT):
        PAGEANT.setObjectName("PAGEANT")
        PAGEANT.setWindowModality(QtCore.Qt.WindowModal)
        PAGEANT.setEnabled(True)
        PAGEANT.resize(600, 480)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred,
                                           QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(60)
        sizePolicy.setVerticalStretch(48)
        sizePolicy.setHeightForWidth(PAGEANT.sizePolicy().hasHeightForWidth())
        PAGEANT.setSizePolicy(sizePolicy)
        PAGEANT.setMinimumSize(QtCore.QSize(600, 480))
        PAGEANT.setMaximumSize(QtCore.QSize(600, 480))

        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(default_font_size)
        font.setBold(True)
        font1 = QtGui.QFont()
        font1.setFamily("Calibri")
        font1.setPointSize(default_font_size - 1)
        font1.setBold(False)
        font2 = QtGui.QFont()
        font2.setFamily("Calibri")
        font2.setPointSize(default_font_size)
        font2.setBold(False)
        font2.setKerning(True)
        font2.setStyleStrategy(QtGui.QFont.PreferDefault)
        font3 = QtGui.QFont()
        font3.setFamily("Calibri")
        font3.setPointSize(default_font_size)
        font3.setBold(False)
        font4 = QtGui.QFont()
        font4.setFamily("Calibri")
        font4.setPointSize(default_font_size)
        font4.setBold(True)
        font4.setKerning(True)
        font4.setStyleStrategy(QtGui.QFont.PreferDefault)

        PAGEANT.setFont(font)
        PAGEANT.setAutoFillBackground(False)
        self.MainWindow = QtWidgets.QWidget(PAGEANT)

        self.MainWindow.setFont(font)
        self.MainWindow.setObjectName("MainWindow")
        self.Function = QtWidgets.QTabWidget(self.MainWindow)
        self.Function.setGeometry(QtCore.QRect(15, 10, 571, 461))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Function.sizePolicy().hasHeightForWidth())
        self.Function.setSizePolicy(sizePolicy)
        self.Function.setMinimumSize(QtCore.QSize(0, 0))

        self.Function.setFont(font)
        self.Function.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.Function.setStyleSheet("")
        self.Function.setIconSize(QtCore.QSize(25, 20))
        self.Function.setDocumentMode(False)
        self.Function.setTabsClosable(False)
        self.Function.setMovable(False)
        self.Function.setTabBarAutoHide(False)
        self.Function.setObjectName("Function")
        self.Basic = QtWidgets.QWidget()

        self.Basic.setFont(font)
        self.Basic.setObjectName("Basic")
        self.i_vcf = QtWidgets.QLineEdit(self.Basic)
        self.i_vcf.setGeometry(QtCore.QRect(79, 141, 365, 25))

        self.i_vcf.setFont(font1)
        self.i_vcf.setObjectName("i_vcf")
        self.b_analyze = QtWidgets.QPushButton(self.Basic)
        self.b_analyze.setGeometry(QtCore.QRect(230, 380, 91, 31))

        self.b_analyze.setFont(font)
        self.b_analyze.setObjectName("b_analyze")
        self.s_vcf = QtWidgets.QToolButton(self.Basic)
        self.s_vcf.setGeometry(QtCore.QRect(440, 140, 41, 27))

        self.s_vcf.setFont(font1)
        self.s_vcf.setObjectName("s_vcf")
        self.l_vcf = QtWidgets.QLabel(self.Basic)
        self.l_vcf.setGeometry(QtCore.QRect(40, 110, 121, 20))

        self.l_vcf.setFont(font)
        self.l_vcf.setObjectName("l_vcf")
        self.l_output = QtWidgets.QLabel(self.Basic)
        self.l_output.setGeometry(QtCore.QRect(40, 290, 171, 20))

        self.l_output.setFont(font)
        self.l_output.setObjectName("l_output")
        self.i_output = QtWidgets.QLineEdit(self.Basic)
        self.i_output.setGeometry(QtCore.QRect(79, 321, 401, 25))

        self.i_output.setFont(font1)
        self.i_output.setObjectName("i_output")
        self.s_output = QtWidgets.QToolButton(self.Basic)
        self.s_output.setGeometry(QtCore.QRect(440, 320, 41, 27))

        self.s_output.setFont(font1)
        self.s_output.setObjectName("s_output")
        self.l_name = QtWidgets.QLabel(self.Basic)
        self.l_name.setGeometry(QtCore.QRect(40, 50, 61, 20))

        self.l_name.setFont(font)
        self.l_name.setObjectName("l_name")
        self.i_name = QtWidgets.QLineEdit(self.Basic)
        self.i_name.setGeometry(QtCore.QRect(100, 50, 121, 25))

        self.i_name.setFont(font1)
        self.i_name.setObjectName("i_name")
        self.l_config = QtWidgets.QLabel(self.Basic)
        self.l_config.setGeometry(QtCore.QRect(40, 200, 201, 25))

        self.l_config.setFont(font)
        self.l_config.setObjectName("l_config")
        self.i_config = QtWidgets.QLineEdit(self.Basic)
        self.i_config.setGeometry(QtCore.QRect(79, 231, 365, 25))

        self.i_config.setFont(font1)
        self.i_config.setObjectName("i_config")
        self.s_config = QtWidgets.QToolButton(self.Basic)
        self.s_config.setGeometry(QtCore.QRect(440, 230, 41, 27))

        self.s_config.setFont(font1)
        self.s_config.setObjectName("s_config")
        self.i_output.raise_()
        self.i_vcf.raise_()
        self.b_analyze.raise_()
        self.l_vcf.raise_()
        self.l_output.raise_()
        self.s_output.raise_()
        self.l_name.raise_()
        self.i_name.raise_()
        self.s_vcf.raise_()
        self.l_config.raise_()
        self.i_config.raise_()
        self.s_config.raise_()
        self.Function.addTab(self.Basic, "")

        self.QC = QtWidgets.QWidget()

        self.QC.setFont(font)
        self.QC.setStyleSheet("")
        self.QC.setObjectName("QC")
        self.c_sample_qc = QtWidgets.QGroupBox(self.QC)
        self.c_sample_qc.setGeometry(QtCore.QRect(20, 20, 521, 351))

        self.c_sample_qc.setFont(font)
        self.c_sample_qc.setCheckable(True)
        self.c_sample_qc.setChecked(True)
        self.c_sample_qc.setObjectName("c_sample_qc")
        self.c_vep = QtWidgets.QCheckBox(self.c_sample_qc)
        self.c_vep.setGeometry(QtCore.QRect(30, 30, 451, 19))

        self.c_vep.setFont(font2)
        self.c_vep.setObjectName("c_vep")
        self.PS = QtWidgets.QGroupBox(self.c_sample_qc)
        self.PS.setEnabled(True)
        self.PS.setGeometry(QtCore.QRect(10, 130, 501, 141))

        self.PS.setFont(font)
        self.PS.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.PS.setObjectName("PS")
        self.l_ps_ref = QtWidgets.QLabel(self.PS)
        self.l_ps_ref.setGeometry(QtCore.QRect(30, 20 + y_offset, 201, 25))

        self.l_ps_ref.setFont(font3)
        self.l_ps_ref.setObjectName("l_ps_ref")
        self.i_ps_ref = QtWidgets.QLineEdit(self.PS)
        self.i_ps_ref.setGeometry(QtCore.QRect(164, 21 + y_offset, 248, 25))

        self.i_ps_ref.setFont(font1)
        self.i_ps_ref.setObjectName("i_ps_ref")
        self.s_ps_ref = QtWidgets.QToolButton(self.PS)
        self.s_ps_ref.setGeometry(QtCore.QRect(410, 20 + y_offset, 41, 27))

        self.s_ps_ref.setFont(font1)
        self.s_ps_ref.setObjectName("s_ps_ref")
        self.l_ps_pop = QtWidgets.QLabel(self.PS)
        self.l_ps_pop.setGeometry(QtCore.QRect(30, 60 + y_offset, 181, 25))

        self.l_ps_pop.setFont(font3)
        self.l_ps_pop.setObjectName("l_ps_pop")
        self.i_pop = QtWidgets.QLineEdit(self.PS)
        self.i_pop.setGeometry(QtCore.QRect(370, 100 + y_offset, 81, 25))

        self.i_pop.setFont(font1)
        self.i_pop.setObjectName("i_pop")
        self.l_pop = QtWidgets.QLabel(self.PS)
        self.l_pop.setGeometry(QtCore.QRect(269, 100 + y_offset, 101, 25))

        self.l_pop.setFont(font3)
        self.l_pop.setObjectName("l_pop")
        self.l_pop_sep = QtWidgets.QLabel(self.PS)
        self.l_pop_sep.setGeometry(QtCore.QRect(30, 100 + y_offset, 61, 25))

        self.l_pop_sep.setFont(font3)
        self.l_pop_sep.setObjectName("l_pop_sep")
        self.s_ps_pop = QtWidgets.QToolButton(self.PS)
        self.s_ps_pop.setGeometry(QtCore.QRect(410, 60 + y_offset, 41, 27))

        self.s_ps_pop.setFont(font1)
        self.s_ps_pop.setObjectName("s_ps_pop")
        self.i_ps_pop = QtWidgets.QLineEdit(self.PS)
        self.i_ps_pop.setGeometry(QtCore.QRect(204, 61 + y_offset, 208, 25))

        self.i_ps_pop.setFont(font1)
        self.i_ps_pop.setObjectName("i_ps_pop")
        self.i_pop_id = QtWidgets.QLineEdit(self.PS)
        self.i_pop_id.setGeometry(QtCore.QRect(210, 100 + y_offset, 41, 25))

        self.i_pop_id.setFont(font1)
        self.i_pop_id.setObjectName("i_pop_id")
        self.l_pop_id = QtWidgets.QLabel(self.PS)
        self.l_pop_id.setGeometry(QtCore.QRect(160, 100 + y_offset, 61, 25))

        self.l_pop_id.setFont(font3)
        self.l_pop_id.setObjectName("l_pop_id")
        self.s_pop_sep = QtWidgets.QComboBox(self.PS)
        self.s_pop_sep.setGeometry(QtCore.QRect(90, 100 + y_offset, 51, 25))

        self.s_pop_sep.setFont(font1)
        self.s_pop_sep.setObjectName("s_pop_sep")
        self.s_pop_sep.addItem("")
        self.s_pop_sep.addItem("")
        self.s_pop_sep.addItem("")
        self.l_ps_ref.raise_()
        self.i_ps_ref.raise_()
        self.s_ps_ref.raise_()
        self.l_ps_pop.raise_()
        self.i_pop.raise_()
        self.l_pop.raise_()
        self.l_pop_sep.raise_()
        self.i_ps_pop.raise_()
        self.i_pop_id.raise_()
        self.l_pop_id.raise_()
        self.s_pop_sep.raise_()
        self.s_ps_pop.raise_()
        self.MAF = QtWidgets.QGroupBox(self.c_sample_qc)
        self.MAF.setEnabled(True)
        self.MAF.setGeometry(QtCore.QRect(10, 70, 501, 61))

        self.MAF.setFont(font)
        self.MAF.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.MAF.setObjectName("MAF")
        self.l_maf_sample = QtWidgets.QLabel(self.MAF)
        self.l_maf_sample.setGeometry(QtCore.QRect(30, 20 + y_offset, 201, 25))

        self.l_maf_sample.setFont(font3)
        self.l_maf_sample.setObjectName("l_maf_sample")
        self.i_maf_sample = QtWidgets.QLineEdit(self.MAF)
        self.i_maf_sample.setGeometry(QtCore.QRect(164, 21 + y_offset, 248, 25))

        self.i_maf_sample.setFont(font1)
        self.i_maf_sample.setObjectName("i_maf_sample")
        self.s_maf_sample = QtWidgets.QToolButton(self.MAF)
        self.s_maf_sample.setGeometry(QtCore.QRect(410, 20 + y_offset, 41, 27))

        self.s_maf_sample.setFont(font1)
        self.s_maf_sample.setObjectName("s_maf_sample")
        self.Concordance = QtWidgets.QGroupBox(self.c_sample_qc)
        self.Concordance.setEnabled(True)
        self.Concordance.setGeometry(QtCore.QRect(10, 280, 501, 61))

        self.Concordance.setFont(font)
        self.Concordance.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.Concordance.setObjectName("Concordance")
        self.l_concord_ref = QtWidgets.QLabel(self.Concordance)
        self.l_concord_ref.setGeometry(QtCore.QRect(30, 20 + y_offset, 201, 25))

        self.l_concord_ref.setFont(font3)
        self.l_concord_ref.setObjectName("l_concord_ref")
        self.i_concord_ref = QtWidgets.QLineEdit(self.Concordance)
        self.i_concord_ref.setGeometry(QtCore.QRect(184, 21 + y_offset, 228, 25))

        self.i_concord_ref.setFont(font1)
        self.i_concord_ref.setObjectName("i_concord_ref")
        self.s_concord_ref = QtWidgets.QToolButton(self.Concordance)
        self.s_concord_ref.setGeometry(QtCore.QRect(410, 20 + y_offset, 41, 27))

        self.s_concord_ref.setFont(font1)
        self.s_concord_ref.setObjectName("s_concord_ref")
        self.b_default_sample_qc = QtWidgets.QPushButton(self.QC)
        self.b_default_sample_qc.setGeometry(QtCore.QRect(230, 380, 91, 31))

        self.b_default_sample_qc.setFont(font)
        self.b_default_sample_qc.setObjectName("b_default_sample_qc")
        self.Function.addTab(self.QC, "")

        self.Qualitative = QtWidgets.QWidget()
        self.Qualitative.setObjectName("Qualitative")
        self.b_default_qual = QtWidgets.QPushButton(self.Qualitative)
        self.b_default_qual.setGeometry(QtCore.QRect(230, 380, 91, 31))

        self.b_default_qual.setFont(font)
        self.b_default_qual.setObjectName("b_default_qual")
        self.g_qual = QtWidgets.QGroupBox(self.Qualitative)
        self.g_qual.setGeometry(QtCore.QRect(20, 20, 521, 351))
        self.g_qual.setObjectName("g_qual")
        self.l_database_qual = QtWidgets.QLabel(self.g_qual)
        self.l_database_qual.setGeometry(QtCore.QRect(30, 80, 201, 25))

        self.l_database_qual.setFont(font)
        self.l_database_qual.setObjectName("l_database_qual")
        self.s_ref_qual = QtWidgets.QToolButton(self.g_qual)
        self.s_ref_qual.setGeometry(QtCore.QRect(430, 230, 41, 27))

        self.s_ref_qual.setFont(font1)
        self.s_ref_qual.setObjectName("s_ref_qual")
        self.i_database_qual = QtWidgets.QLineEdit(self.g_qual)
        self.i_database_qual.setGeometry(QtCore.QRect(69, 111, 365, 25))

        self.i_database_qual.setFont(font1)
        self.i_database_qual.setObjectName("i_database_qual")
        self.i_ref_qual = QtWidgets.QLineEdit(self.g_qual)
        self.i_ref_qual.setGeometry(QtCore.QRect(69, 231, 365, 25))

        self.i_ref_qual.setFont(font1)
        self.i_ref_qual.setObjectName("i_ref_qual")
        self.s_database_qual = QtWidgets.QToolButton(self.g_qual)
        self.s_database_qual.setGeometry(QtCore.QRect(430, 110, 41, 27))

        self.s_database_qual.setFont(font1)
        self.s_database_qual.setObjectName("s_database_qual")
        self.l_ref_qual = QtWidgets.QLabel(self.g_qual)
        self.l_ref_qual.setGeometry(QtCore.QRect(30, 200, 201, 25))

        self.l_ref_qual.setFont(font)
        self.l_ref_qual.setObjectName("l_ref_qual")
        self.g_qual.raise_()
        self.b_default_qual.raise_()
        self.Function.addTab(self.Qualitative, "")

        self.Quantitative = QtWidgets.QWidget()

        self.Quantitative.setFont(font)
        self.Quantitative.setStyleSheet("")
        self.Quantitative.setObjectName("Quantitative")
        self.b_default_quan = QtWidgets.QPushButton(self.Quantitative)
        self.b_default_quan.setGeometry(QtCore.QRect(230, 380, 91, 31))

        self.b_default_quan.setFont(font)
        self.b_default_quan.setObjectName("b_default_quan")
        self.g_quan = QtWidgets.QGroupBox(self.Quantitative)
        self.g_quan.setGeometry(QtCore.QRect(20, 20, 521, 351))
        self.g_quan.setObjectName("g_quan")
        self.l_database_quan = QtWidgets.QLabel(self.g_quan)
        self.l_database_quan.setGeometry(QtCore.QRect(30, 50, 201, 25))

        self.l_database_quan.setFont(font)
        self.l_database_quan.setObjectName("l_database_quan")
        self.s_ref_quan = QtWidgets.QToolButton(self.g_quan)
        self.s_ref_quan.setGeometry(QtCore.QRect(430, 170, 41, 27))

        self.s_ref_quan.setFont(font1)
        self.s_ref_quan.setObjectName("s_ref_quan")
        self.i_database_quan = QtWidgets.QLineEdit(self.g_quan)
        self.i_database_quan.setGeometry(QtCore.QRect(69, 81, 365, 25))

        self.i_database_quan.setFont(font1)
        self.i_database_quan.setObjectName("i_database_quan")
        self.i_ref_quan = QtWidgets.QLineEdit(self.g_quan)
        self.i_ref_quan.setGeometry(QtCore.QRect(69, 171, 365, 25))

        self.i_ref_quan.setFont(font1)
        self.i_ref_quan.setObjectName("i_ref_quan")
        self.s_database_quan = QtWidgets.QToolButton(self.g_quan)
        self.s_database_quan.setGeometry(QtCore.QRect(430, 80, 41, 27))

        self.s_database_quan.setFont(font1)
        self.s_database_quan.setObjectName("s_database_quan")
        self.l_ref_quan = QtWidgets.QLabel(self.g_quan)
        self.l_ref_quan.setGeometry(QtCore.QRect(30, 140, 201, 25))

        self.l_ref_quan.setFont(font)
        self.l_ref_quan.setObjectName("l_ref_quan")
        self.l_quan_pop = QtWidgets.QLabel(self.g_quan)
        self.l_quan_pop.setGeometry(QtCore.QRect(30, 230, 371, 25))

        self.l_quan_pop.setFont(font)
        self.l_quan_pop.setObjectName("l_quan_pop")
        self.s_quan_pop = QtWidgets.QComboBox(self.g_quan)
        self.s_quan_pop.setGeometry(QtCore.QRect(70, 260, 401, 25))

        self.s_quan_pop.setFont(font1)
        self.s_quan_pop.setObjectName("s_quan_pop")
        self.s_quan_pop.addItem("")
        self.s_quan_pop.addItem("")
        self.s_quan_pop.addItem("")
        self.s_quan_pop.addItem("")
        self.s_quan_pop.addItem("")
        self.s_quan_pop.addItem("")
        self.Function.addTab(self.Quantitative, "")

        self.Query_database = QtWidgets.QWidget()

        self.Query_database.setFont(font)
        self.Query_database.setObjectName("Query_database")
        self.c_query_database = QtWidgets.QGroupBox(self.Query_database)
        self.c_query_database.setGeometry(QtCore.QRect(30, 20, 501, 351))

        self.c_query_database.setFont(font)
        self.c_query_database.setCheckable(True)
        self.c_query_database.setObjectName("c_query_database")
        self.c_clinvar = QtWidgets.QCheckBox(self.c_query_database)
        self.c_clinvar.setGeometry(QtCore.QRect(40, 60, 451, 19))

        self.c_clinvar.setFont(font4)
        self.c_clinvar.setChecked(True)
        self.c_clinvar.setObjectName("c_clinvar")
        self.c_pharmgkb = QtWidgets.QCheckBox(self.c_query_database)
        self.c_pharmgkb.setGeometry(QtCore.QRect(40, 110, 451, 19))

        self.c_pharmgkb.setFont(font4)
        self.c_pharmgkb.setChecked(True)
        self.c_pharmgkb.setObjectName("c_pharmgkb")
        self.c_otherdb = QtWidgets.QGroupBox(self.c_query_database)
        self.c_otherdb.setGeometry(QtCore.QRect(32, 160, 441, 141))
        self.c_otherdb.setCheckable(True)
        self.c_otherdb.setChecked(True)
        self.c_otherdb.setObjectName("c_otherdb")
        self.l_otherdb = QtWidgets.QLabel(self.c_otherdb)
        self.l_otherdb.setGeometry(QtCore.QRect(27, 40, 201, 25))

        self.l_otherdb.setFont(font3)
        self.l_otherdb.setObjectName("l_otherdb")
        self.i_otherdb = QtWidgets.QLineEdit(self.c_otherdb)
        self.i_otherdb.setGeometry(QtCore.QRect(31, 81, 341, 25))

        self.i_otherdb.setFont(font1)
        self.i_otherdb.setObjectName("i_otherdb")
        self.s_otherdb = QtWidgets.QToolButton(self.c_otherdb)
        self.s_otherdb.setGeometry(QtCore.QRect(370, 80, 41, 27))

        self.s_otherdb.setFont(font1)
        self.s_otherdb.setObjectName("s_otherdb")
        self.b_default_query_db = QtWidgets.QPushButton(self.Query_database)
        self.b_default_query_db.setGeometry(QtCore.QRect(230, 380, 91, 31))

        self.b_default_query_db.setFont(font)
        self.b_default_query_db.setObjectName("b_default_query_db")
        self.Function.addTab(self.Query_database, "")
        self.QR_code = QtWidgets.QWidget()

        self.QR_code.setFont(font)
        self.QR_code.setObjectName("QR_code")
        self.c_qr_code = QtWidgets.QGroupBox(self.QR_code)
        self.c_qr_code.setGeometry(QtCore.QRect(30, 20, 501, 351))

        self.c_qr_code.setFont(font)
        self.c_qr_code.setCheckable(True)
        self.c_qr_code.setObjectName("c_qr_code")

        self.l_qr_give = QtWidgets.QLabel(self.c_qr_code)
        self.l_qr_give.setGeometry(QtCore.QRect(30, 180, 281, 25))

        self.l_qr_give.setFont(font)
        self.l_qr_give.setObjectName("l_qr_give")
        self.i_qr_give = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_qr_give.setGeometry(QtCore.QRect(50, 211, 365, 25))

        self.i_qr_give.setFont(font1)
        self.i_qr_give.setObjectName("i_qr_give")
        self.s_qr_give = QtWidgets.QToolButton(self.c_qr_code)
        self.s_qr_give.setGeometry(QtCore.QRect(410, 210, 41, 27))

        self.s_qr_give.setFont(font1)
        self.s_qr_give.setObjectName("s_qr_give")
        self.l_qr_dir = QtWidgets.QLabel(self.c_qr_code)
        self.l_qr_dir.setGeometry(QtCore.QRect(30, 40, 281, 25))

        self.l_qr_dir.setFont(font)
        self.l_qr_dir.setObjectName("l_qr_dir")
        self.i_qr_dir = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_qr_dir.setGeometry(QtCore.QRect(50, 71, 365, 25))

        self.i_qr_dir.setFont(font1)
        self.i_qr_dir.setObjectName("i_qr_dir")
        self.s_qr_dir = QtWidgets.QToolButton(self.c_qr_code)
        self.s_qr_dir.setGeometry(QtCore.QRect(410, 70, 41, 27))

        self.s_qr_dir.setFont(font1)
        self.s_qr_dir.setObjectName("s_qr_dir")
        self.l_qr_user = QtWidgets.QLabel(self.c_qr_code)
        self.l_qr_user.setGeometry(QtCore.QRect(30, 250, 251, 25))

        self.l_qr_user.setFont(font)
        self.l_qr_user.setObjectName("l_qr_user")
        self.i_qr_user = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_qr_user.setGeometry(QtCore.QRect(50, 281, 365, 25))

        self.i_qr_user.setFont(font1)
        self.i_qr_user.setObjectName("i_qr_user")
        self.s_qr_user = QtWidgets.QToolButton(self.c_qr_code)
        self.s_qr_user.setGeometry(QtCore.QRect(410, 280, 41, 27))

        self.s_qr_user.setFont(font1)
        self.s_qr_user.setObjectName("s_qr_user")
        self.l_qr_snps = QtWidgets.QLabel(self.c_qr_code)
        self.l_qr_snps.setGeometry(QtCore.QRect(30, 110, 281, 25))

        self.l_qr_snps.setFont(font)
        self.l_qr_snps.setObjectName("l_qr_snps")
        self.i_qr_snps = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_qr_snps.setGeometry(QtCore.QRect(50, 141, 365, 25))

        self.i_qr_snps.setFont(font1)
        self.i_qr_snps.setObjectName("i_qr_snps")
        self.s_qr_snps = QtWidgets.QToolButton(self.c_qr_code)
        self.s_qr_snps.setGeometry(QtCore.QRect(410, 140, 41, 27))

        self.s_qr_snps.setFont(font1)
        self.s_qr_snps.setObjectName("s_qr_snps")

        self.b_default_qr_code = QtWidgets.QPushButton(self.QR_code)
        self.b_default_qr_code.setGeometry(QtCore.QRect(230, 380, 91, 31))

        self.b_default_qr_code.setFont(font)
        self.b_default_qr_code.setObjectName("b_default_qr_code")
        self.Function.addTab(self.QR_code, "")

        self.API = QtWidgets.QWidget()
        self.API.setObjectName("API")
        self.UMAP = QtWidgets.QGroupBox(self.API)
        self.UMAP.setEnabled(True)
        self.UMAP.setGeometry(QtCore.QRect(30, 10, 501, 161))

        self.UMAP.setFont(font)
        self.UMAP.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.UMAP.setObjectName("UMAP")
        self.l_umap_ref = QtWidgets.QLabel(self.UMAP)
        self.l_umap_ref.setGeometry(QtCore.QRect(30, 31, 131, 25))

        self.l_umap_ref.setFont(font2)
        self.l_umap_ref.setObjectName("l_umap_ref")
        self.i_umap_ref = QtWidgets.QLineEdit(self.UMAP)
        self.i_umap_ref.setGeometry(QtCore.QRect(164, 32, 248, 25))

        self.i_umap_ref.setFont(font1)
        self.i_umap_ref.setObjectName("i_umap_ref")
        self.s_umap_ref = QtWidgets.QToolButton(self.UMAP)
        self.s_umap_ref.setGeometry(QtCore.QRect(410, 31, 41, 27))

        self.s_umap_ref.setFont(font1)
        self.s_umap_ref.setObjectName("s_umap_ref")
        self.l_umap_sample = QtWidgets.QLabel(self.UMAP)
        self.l_umap_sample.setGeometry(QtCore.QRect(30, 61, 131, 25))

        self.l_umap_sample.setFont(font2)
        self.l_umap_sample.setObjectName("l_umap_sample")
        self.i_umap_sample = QtWidgets.QLineEdit(self.UMAP)
        self.i_umap_sample.setGeometry(QtCore.QRect(164, 62, 248, 25))

        self.i_umap_sample.setFont(font1)
        self.i_umap_sample.setObjectName("i_umap_sample")
        self.s_umap_sample = QtWidgets.QToolButton(self.UMAP)
        self.s_umap_sample.setGeometry(QtCore.QRect(410, 61, 41, 27))

        self.s_umap_sample.setFont(font1)
        self.s_umap_sample.setObjectName("s_umap_sample")
        self.l_umap_metadata = QtWidgets.QLabel(self.UMAP)
        self.l_umap_metadata.setGeometry(QtCore.QRect(30, 90, 131, 25))

        self.l_umap_metadata.setFont(font2)
        self.l_umap_metadata.setObjectName("l_umap_metadata")
        self.s_umap_metadata = QtWidgets.QToolButton(self.UMAP)
        self.s_umap_metadata.setGeometry(QtCore.QRect(410, 90, 41, 27))

        self.s_umap_metadata.setFont(font1)
        self.s_umap_metadata.setObjectName("s_umap_metadata")
        self.i_umap_metadata = QtWidgets.QLineEdit(self.UMAP)
        self.i_umap_metadata.setGeometry(QtCore.QRect(164, 91, 248, 25))

        self.i_umap_metadata.setFont(font1)
        self.i_umap_metadata.setObjectName("i_umap_metadata")
        self.b_run_umap = QtWidgets.QPushButton(self.UMAP)
        self.b_run_umap.setGeometry(QtCore.QRect(170, 120, 181, 31))

        self.b_run_umap.setFont(font)
        self.b_run_umap.setObjectName("b_run_umap")
        self.add_rsid = QtWidgets.QGroupBox(self.API)
        self.add_rsid.setEnabled(True)
        self.add_rsid.setGeometry(QtCore.QRect(30, 180, 501, 111))

        self.add_rsid.setFont(font)
        self.add_rsid.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.add_rsid.setObjectName("add_rsid")
        self.l_add_rsid_file = QtWidgets.QLabel(self.add_rsid)
        self.l_add_rsid_file.setGeometry(QtCore.QRect(30, 31, 131, 25))

        self.l_add_rsid_file.setFont(font2)
        self.l_add_rsid_file.setObjectName("l_add_rsid_file")
        self.i_add_rsid_file = QtWidgets.QLineEdit(self.add_rsid)
        self.i_add_rsid_file.setGeometry(QtCore.QRect(121, 32, 221, 25))

        self.i_add_rsid_file.setFont(font1)
        self.i_add_rsid_file.setObjectName("i_add_rsid_file")
        self.s_add_rsid_file = QtWidgets.QToolButton(self.add_rsid)
        self.s_add_rsid_file.setGeometry(QtCore.QRect(340, 31, 41, 27))

        self.s_add_rsid_file.setFont(font1)
        self.s_add_rsid_file.setObjectName("s_add_rsid_file")
        self.l_add_rsid_db = QtWidgets.QLabel(self.add_rsid)
        self.l_add_rsid_db.setGeometry(QtCore.QRect(30, 70, 131, 25))

        self.l_add_rsid_db.setFont(font2)
        self.l_add_rsid_db.setObjectName("l_add_rsid_db")
        self.i_add_rsid_db = QtWidgets.QLineEdit(self.add_rsid)
        self.i_add_rsid_db.setGeometry(QtCore.QRect(121, 71, 221, 25))
        self.s_add_rsid_db = QtWidgets.QToolButton(self.add_rsid)
        self.s_add_rsid_db.setGeometry(QtCore.QRect(340, 70, 41, 27))

        self.s_add_rsid_db.setFont(font1)
        self.s_add_rsid_db.setObjectName("s_add_rsid_db")

        self.i_add_rsid_db.setFont(font1)
        self.i_add_rsid_db.setObjectName("i_add_rsid_db")
        self.b_run_add_rsid = QtWidgets.QPushButton(self.add_rsid)
        self.b_run_add_rsid.setGeometry(QtCore.QRect(390, 33, 101, 61))

        self.b_run_add_rsid.setFont(font)
        self.b_run_add_rsid.setObjectName("b_run_add_rsid")
        self.qr_code = QtWidgets.QGroupBox(self.API)
        self.qr_code.setEnabled(True)
        self.qr_code.setGeometry(QtCore.QRect(30, 300, 501, 111))

        self.qr_code.setFont(font)
        self.qr_code.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.qr_code.setObjectName("qr_code")
        self.l_qr_code_snp = QtWidgets.QLabel(self.qr_code)
        self.l_qr_code_snp.setGeometry(QtCore.QRect(30, 31, 121, 25))

        self.l_qr_code_snp.setFont(font2)
        self.l_qr_code_snp.setObjectName("l_qr_code_snp")
        self.i_qr_code_snp = QtWidgets.QLineEdit(self.qr_code)
        self.i_qr_code_snp.setGeometry(QtCore.QRect(121, 32, 221, 25))

        self.i_qr_code_snp.setFont(font1)
        self.i_qr_code_snp.setObjectName("i_qr_code_snp")
        self.s_qr_code_snp = QtWidgets.QToolButton(self.qr_code)
        self.s_qr_code_snp.setGeometry(QtCore.QRect(340, 31, 41, 27))

        self.s_qr_code_snp.setFont(font1)
        self.s_qr_code_snp.setObjectName("s_qr_code_snp")
        self.b_run_qr_code = QtWidgets.QPushButton(self.qr_code)
        self.b_run_qr_code.setGeometry(QtCore.QRect(390, 33, 101, 61))

        self.b_run_qr_code.setFont(font)
        self.b_run_qr_code.setObjectName("b_run_qr_code")
        self.s_qr_code_key = QtWidgets.QToolButton(self.qr_code)
        self.s_qr_code_key.setGeometry(QtCore.QRect(340, 70, 41, 27))

        self.s_qr_code_key.setFont(font1)
        self.s_qr_code_key.setObjectName("s_qr_code_key")
        self.l_qr_code_key = QtWidgets.QLabel(self.qr_code)
        self.l_qr_code_key.setGeometry(QtCore.QRect(30, 70, 121, 25))

        self.l_qr_code_key.setFont(font2)
        self.l_qr_code_key.setObjectName("l_qr_code_key")
        self.i_qr_code_key = QtWidgets.QLineEdit(self.qr_code)
        self.i_qr_code_key.setGeometry(QtCore.QRect(121, 71, 221, 25))

        self.i_qr_code_key.setFont(font1)
        self.i_qr_code_key.setObjectName("i_qr_code_key")
        self.Function.addTab(self.API, "")

        PAGEANT.setCentralWidget(self.MainWindow)

        self.retranslateUi(PAGEANT)
        self.Function.setCurrentIndex(0)
        self.s_pop_sep.setCurrentIndex(0)
        self.s_quan_pop.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(PAGEANT)

        self.l_database_quan.raise_()
        self.i_database_quan.raise_()
        self.i_ref_quan.raise_()
        self.s_database_quan.raise_()
        self.l_ref_quan.raise_()
        self.l_quan_pop.raise_()
        self.s_ref_quan.raise_()
        self.s_quan_pop.raise_()
        self.l_database_qual.raise_()
        self.i_database_qual.raise_()
        self.i_ref_qual.raise_()
        self.s_database_qual.raise_()
        self.l_ref_qual.raise_()
        self.s_ref_qual.raise_()
        self.g_qual.raise_()
        self.b_default_qual.raise_()
        self.i_output.raise_()
        self.i_vcf.raise_()
        self.b_analyze.raise_()
        self.l_vcf.raise_()
        self.l_output.raise_()
        self.s_output.raise_()
        self.l_name.raise_()
        self.i_name.raise_()
        self.s_vcf.raise_()
        self.l_config.raise_()
        self.i_config.raise_()
        self.s_config.raise_()
        self.l_ps_ref.raise_()
        self.i_ps_ref.raise_()
        self.s_ps_ref.raise_()
        self.l_ps_pop.raise_()
        self.i_pop.raise_()
        self.l_pop.raise_()
        self.l_pop_sep.raise_()
        self.i_ps_pop.raise_()
        self.i_pop_id.raise_()
        self.l_pop_id.raise_()
        self.s_pop_sep.raise_()
        self.s_ps_pop.raise_()

    def retranslateUi(self, PAGEANT):
        _translate = QtCore.QCoreApplication.translate
        PAGEANT.setWindowTitle(_translate("PAGEANT", f"PAGEANT ({main.version})"))
        self.i_vcf.setToolTip(_translate("PAGEANT", "Select the genotype file that you want to analyze"))
        self.i_vcf.setText(_translate("PAGEANT", "./personal_genome/HG001.vcf.gz"))
        self.i_vcf.setPlaceholderText(_translate("PAGEANT", "Select the genotype file"))
        self.b_analyze.setToolTip(_translate("PAGEANT", "Start analyze"))
        self.b_analyze.setText(_translate("PAGEANT", "Analyze"))
        self.s_vcf.setText(_translate("PAGEANT", "..."))
        self.l_vcf.setToolTip(_translate("PAGEANT", "Input the genotype file that you want to analyze"))
        self.l_vcf.setText(_translate("PAGEANT", "Genotype file"))
        self.l_output.setText(_translate("PAGEANT", "Ouput Directory"))
        self.i_output.setToolTip(_translate("PAGEANT", "Specify the output directory"))
        self.i_output.setText(_translate("PAGEANT", os.path.join(os.getcwd(), 'output')))
        self.i_output.setPlaceholderText(_translate("PAGEANT", 'Specify the output directory'))
        self.s_output.setText(_translate("PAGEANT", "..."))
        self.l_name.setText(_translate("PAGEANT", "Name"))
        self.i_name.setToolTip(_translate("PAGEANT", "Input your name here"))
        self.i_name.setText(_translate("PAGEANT", "HG001"))
        self.l_config.setText(_translate("PAGEANT", "Configuration file"))
        self.i_config.setToolTip(_translate("PAGEANT", "Specify the configuration file"))
        self.i_config.setText(_translate("PAGEANT", "./bin/config.ini"))
        self.i_config.setPlaceholderText(_translate("PAGEANT", "Specify the configuration file"))
        self.s_config.setText(_translate("PAGEANT", "..."))
        self.Function.setTabText(self.Function.indexOf(self.Basic), _translate("PAGEANT", "I/O"))

        self.b_default_qual.setToolTip(_translate("PAGEANT", "Set default values"))
        self.b_default_qual.setText(_translate("PAGEANT", "Default"))
        self.g_qual.setTitle(_translate("PAGEANT", "Qualitative traits"))
        self.l_database_qual.setText(_translate("PAGEANT", "Database directory"))
        self.s_ref_qual.setText(_translate("PAGEANT", "..."))
        self.i_database_qual.setToolTip(_translate("PAGEANT", "Specify directory of database for qualitative traits"))
        self.i_database_qual.setText(_translate("PAGEANT", "./algorithm/Qualitative"))
        self.i_database_qual.setPlaceholderText(_translate("PAGEANT", "Specify database for qualitative traits"))
        self.i_ref_qual.setToolTip(_translate("PAGEANT", "Specify directory of reference for qualitative traits"))
        self.i_ref_qual.setText(_translate("PAGEANT", "./population_genome"))
        self.i_ref_qual.setPlaceholderText(_translate("PAGEANT", "Specify reference for qualitative traits"))
        self.s_database_qual.setText(_translate("PAGEANT", "..."))
        self.l_ref_qual.setText(_translate("PAGEANT", "Reference population directory"))
        self.Function.setTabText(self.Function.indexOf(self.Qualitative), _translate("PAGEANT", "Qualitative"))
        self.b_default_quan.setToolTip(_translate("PAGEANT", "Set default values"))
        self.b_default_quan.setText(_translate("PAGEANT", "Default"))
        self.g_quan.setTitle(_translate("PAGEANT", "Quantitative traits"))
        self.l_database_quan.setText(_translate("PAGEANT", "Database directory"))
        self.s_ref_quan.setText(_translate("PAGEANT", "..."))
        self.i_database_quan.setToolTip(_translate("PAGEANT", "Specify directory of database for quantitative traits"))
        self.i_database_quan.setText(_translate("PAGEANT", "./algorithm/Quantitative"))
        self.i_database_quan.setPlaceholderText(_translate("PAGEANT", "Specify database for quantitative traits"))
        self.i_ref_quan.setToolTip(_translate("PAGEANT", "Specify directory of reference for quantitative traits"))
        self.i_ref_quan.setText(_translate("PAGEANT", "./population_genome"))
        self.i_ref_quan.setPlaceholderText(_translate("PAGEANT", "Specify reference for quantitative traits"))
        self.s_database_quan.setText(_translate("PAGEANT", "..."))
        self.l_ref_quan.setText(_translate("PAGEANT", "Reference population directory"))
        self.l_quan_pop.setText(_translate("PAGEANT", "Reference population ethnical group"))
        self.s_quan_pop.setPlaceholderText(_translate("PAGEANT", "Select the separator of population data"))
        self.s_quan_pop.setItemText(0, _translate("PAGEANT", "All"))
        self.s_quan_pop.setItemText(1, _translate("PAGEANT", "EUR"))
        self.s_quan_pop.setItemText(2, _translate("PAGEANT", "EAS"))
        self.s_quan_pop.setItemText(3, _translate("PAGEANT", "AMR"))
        self.s_quan_pop.setItemText(4, _translate("PAGEANT", "SAS"))
        self.s_quan_pop.setItemText(5, _translate("PAGEANT", "AFR"))
        self.Function.setTabText(self.Function.indexOf(self.Quantitative), _translate("PAGEANT", "Quantitative"))

        # self.c_sample_qc.setToolTip(_translate("PAGEANT", "Activate function of sample QC"))
        self.c_sample_qc.setTitle(_translate("PAGEANT", "Sample QC"))
        self.c_vep.setToolTip(_translate("PAGEANT", "Using VEP to inspect the variants\' type in sample data"))
        self.c_vep.setText(_translate("PAGEANT", "Use VEP (Slow with too many variants or bad network!)"))
        self.PS.setTitle(_translate("PAGEANT", "Population stratification analysis"))
        self.l_ps_ref.setText(_translate("PAGEANT", "Reference data"))
        self.i_ps_ref.setToolTip(_translate("PAGEANT", "Specify the reference data for ps"))
        self.i_ps_ref.setText(_translate("PAGEANT", "./personal_genome/g1k.vcf.gz"))
        self.i_ps_ref.setPlaceholderText(_translate("PAGEANT", "Select principal component analysis reference data"))
        self.s_ps_ref.setText(_translate("PAGEANT", "..."))
        self.l_ps_pop.setText(_translate("PAGEANT", "Reference population data"))
        self.i_pop.setToolTip(_translate("PAGEANT", "Separator in GWAS data"))
        self.i_pop.setText(_translate("PAGEANT", "population"))
        self.i_pop.setPlaceholderText(_translate("PAGEANT", "Input the column\' name of population"))
        self.l_pop.setText(_translate("PAGEANT", "Population Col."))
        self.l_pop_sep.setText(_translate("PAGEANT", "Text Sep."))
        self.s_ps_pop.setText(_translate("PAGEANT", "..."))
        self.i_ps_pop.setToolTip(_translate("PAGEANT", "Specify the reference population data for "
                                                       "population stratification analysis"))
        self.i_ps_pop.setText(_translate("PAGEANT", "./personal_genome/g1k_samples.txt"))
        self.i_ps_pop.setPlaceholderText(_translate("PAGEANT", "Select principal component analysis population data"))
        self.i_pop_id.setToolTip(_translate("PAGEANT", "Separator in *.snps.ref"))
        self.i_pop_id.setText(_translate("PAGEANT", "IID"))
        self.i_pop_id.setPlaceholderText(_translate("PAGEANT", "Input the column\' name of ID"))
        self.l_pop_id.setText(_translate("PAGEANT", "ID Col."))
        self.s_pop_sep.setPlaceholderText(_translate("PAGEANT", "Select the separator of population data"))
        self.s_pop_sep.setItemText(0, _translate("PAGEANT", "\\t"))
        self.s_pop_sep.setItemText(1, _translate("PAGEANT", ","))
        self.s_pop_sep.setItemText(2, _translate("PAGEANT", "SPACE"))
        self.MAF.setTitle(_translate("PAGEANT", "MAF plot"))
        self.l_maf_sample.setText(_translate("PAGEANT", "MAF reference data"))
        self.i_maf_sample.setToolTip(_translate("PAGEANT", "Specify the reference data for MAF QC"))
        self.i_maf_sample.setText(_translate("PAGEANT", "./personal_genome/g1k.vcf.gz"))
        self.i_maf_sample.setPlaceholderText(_translate("PAGEANT", "Select MAF reference data"))
        self.s_maf_sample.setText(_translate("PAGEANT", "..."))
        self.Concordance.setTitle(_translate("PAGEANT", "Concordance"))
        self.l_concord_ref.setText(_translate("PAGEANT", "Concordance ref. data"))
        self.i_concord_ref.setToolTip(_translate("PAGEANT", "Specify the reference data for concordance check"))
        self.i_concord_ref.setText(_translate("PAGEANT", "./personal_genome/concordance.vcf.gz"))
        self.i_concord_ref.setPlaceholderText(_translate("PAGEANT", "Select concordance analysis reference data"))
        self.s_concord_ref.setText(_translate("PAGEANT", "..."))
        self.b_default_sample_qc.setToolTip(_translate("PAGEANT", "Set default values"))
        self.b_default_sample_qc.setText(_translate("PAGEANT", "Default"))
        self.Function.setTabText(self.Function.indexOf(self.QC), _translate("PAGEANT", "QA/QC"))
        # self.c_query_database.setToolTip(_translate("PAGEANT", "Activate function of querying database"))
        self.c_query_database.setTitle(_translate("PAGEANT", "Query database"))
        self.c_clinvar.setToolTip(_translate("PAGEANT", "Query Clinvar Database"))
        self.c_clinvar.setText(_translate("PAGEANT", "Query Clinvar Database"))
        self.c_pharmgkb.setToolTip(_translate("PAGEANT", "Query PharmGKB Database"))
        self.c_pharmgkb.setText(_translate("PAGEANT", "Query PharmGKB Database"))
        self.b_default_query_db.setToolTip(_translate("PAGEANT", "Set default values"))
        self.b_default_query_db.setText(_translate("PAGEANT", "Default"))
        self.c_otherdb.setTitle(_translate("PAGEANT", "Query third-party database"))
        self.l_otherdb.setText(_translate("PAGEANT", "Third-party database"))
        self.i_otherdb.setToolTip(_translate("PAGEANT", "Specify the third-party database file"))
        self.i_otherdb.setText(_translate("PAGEANT", "./algorithm/Query_db/Phewas_Catalog_part.tsv"))
        self.i_otherdb.setPlaceholderText(_translate("PAGEANT", "Select third-party database file"))
        self.s_otherdb.setText(_translate("PAGEANT", "..."))
        self.Function.setTabText(self.Function.indexOf(self.Query_database), _translate("PAGEANT", "Query DB"))

        # self.c_qr_code.setToolTip(_translate("PAGEANT", "Activate function of generating QR code"))
        self.c_qr_code.setTitle(_translate("PAGEANT", "Generate QR code"))

        self.l_qr_give.setText(_translate("PAGEANT", "SNP list QR code (for user to scan)"))
        self.i_qr_give.setToolTip(_translate("PAGEANT", "Specify the path for saving SNP list QR code"))
        self.i_qr_give.setText(_translate("PAGEANT", "./output/qr_code/snp_qr_code.png"))
        self.s_qr_give.setText(_translate("PAGEANT", "..."))
        self.l_qr_dir.setText(_translate("PAGEANT", "Directory"))
        self.i_qr_dir.setToolTip(_translate("PAGEANT", "Specify the directory where needed files can be found"))
        self.i_qr_dir.setText(_translate("PAGEANT", "./output/qr_code"))
        self.s_qr_dir.setText(_translate("PAGEANT", "..."))
        self.l_qr_user.setText(_translate("PAGEANT", "User QR code (with encrypted genotype)"))
        self.i_qr_user.setToolTip(_translate("PAGEANT", "Specify the path for saving User QR code"))
        self.i_qr_user.setText(_translate("PAGEANT", "./output/qr_code/user_qr_code.png"))
        self.s_qr_user.setText(_translate("PAGEANT", "..."))
        self.l_qr_snps.setText(_translate("PAGEANT", "SNP list file"))
        self.i_qr_snps.setToolTip(_translate("PAGEANT", "Specify the text which include the needed SNPs list"))
        self.i_qr_snps.setText(_translate("PAGEANT", "./personal_genome/fingerprint_snps.txt"))
        self.s_qr_snps.setText(_translate("PAGEANT", "..."))

        self.b_default_qr_code.setToolTip(_translate("PAGEANT", "Set default values"))
        self.b_default_qr_code.setText(_translate("PAGEANT", "Default"))
        self.Function.setTabText(self.Function.indexOf(self.QR_code), _translate("PAGEANT", "QR code"))

        self.UMAP.setTitle(_translate("PAGEANT", "The API to generate UMAP and PCA plot"))
        self.l_umap_ref.setText(_translate("PAGEANT", "Reference genome"))
        self.i_umap_ref.setToolTip(_translate("PAGEANT", "Select the reference file for PCA & UMAP plot"))
        self.i_umap_ref.setText(_translate("PAGEANT", "./personal_genome/g1k.vcf.gz"))
        self.i_umap_ref.setPlaceholderText(_translate("PAGEANT", "Select reference genome"))
        self.s_umap_ref.setText(_translate("PAGEANT", "..."))
        self.l_umap_sample.setText(_translate("PAGEANT", "Sample genome"))
        self.i_umap_sample.setToolTip(_translate("PAGEANT", "Select the sample file"))
        self.i_umap_sample.setText(_translate("PAGEANT", "./personal_genome/HG001.vcf.gz"))
        self.i_umap_sample.setPlaceholderText(_translate("PAGEANT", "Select sample genome"))
        self.s_umap_sample.setText(_translate("PAGEANT", "..."))
        self.l_umap_metadata.setText(_translate("PAGEANT", "Metadata"))
        self.s_umap_metadata.setText(_translate("PAGEANT", "..."))
        self.i_umap_metadata.setToolTip(_translate("PAGEANT", "Select the metadata for reference genotype data"))
        self.i_umap_metadata.setText(_translate("PAGEANT", "./personal_genome/g1k_samples.txt"))
        self.i_umap_metadata.setPlaceholderText(_translate("PAGEANT", "Select metadata"))
        self.b_run_umap.setToolTip(_translate("PAGEANT", "RUN PCA and UMAP plot"))
        self.b_run_umap.setText(_translate("PAGEANT", "RUN PCA and UMAP"))
        self.add_rsid.setTitle(_translate("PAGEANT", "The API for add rsID to GWAS"))
        self.l_add_rsid_file.setText(_translate("PAGEANT", "GWAS file"))
        self.i_add_rsid_file.setToolTip(_translate("PAGEANT", "Select the GWAS file"))
        self.i_add_rsid_file.setText(_translate("PAGEANT", "./add_rsid/test.tsv"))
        self.i_add_rsid_file.setPlaceholderText(_translate("PAGEANT", "Specific the format"))
        self.s_add_rsid_file.setText(_translate("PAGEANT", "..."))
        self.l_add_rsid_db.setText(_translate("PAGEANT", "dbSNP file"))
        self.i_add_rsid_db.setToolTip(_translate("PAGEANT", "Specify the path of dbSNP file"))
        self.i_add_rsid_db.setText(_translate("PAGEANT", "./add_rsid/rsids-v150-hg19.tsv.gz"))
        self.i_add_rsid_db.setPlaceholderText(_translate("PAGEANT", "Specify the path of dbSNP file"))
        self.s_add_rsid_db.setText(_translate("PAGEANT", "..."))
        self.b_run_add_rsid.setToolTip(_translate("PAGEANT", "Start analyze"))
        self.b_run_add_rsid.setText(_translate("PAGEANT", "Add rsID \nto GWAS file"))
        self.qr_code.setTitle(_translate("PAGEANT", "The API for generating SNP QR code"))
        self.l_qr_code_snp.setText(_translate("PAGEANT", "SNP list"))
        self.i_qr_code_snp.setToolTip(_translate("PAGEANT", "Select the text which include the needed SNPs list"))
        self.i_qr_code_snp.setText(_translate("PAGEANT", "./personal_genome/fingerprint_snps.txt"))
        self.i_qr_code_snp.setPlaceholderText(_translate("PAGEANT", "Select the file which contains public key"))
        self.s_qr_code_snp.setText(_translate("PAGEANT", "..."))
        self.b_run_qr_code.setToolTip(_translate("PAGEANT", "Generate SNP QR CODE"))
        self.b_run_qr_code.setText(_translate("PAGEANT", "Create \nQR Code"))
        self.s_qr_code_key.setText(_translate("PAGEANT", "..."))
        self.l_qr_code_key.setText(_translate("PAGEANT", "Public key"))
        self.i_qr_code_key.setToolTip(_translate("PAGEANT", "Specify directory of database"))
        self.i_qr_code_key.setText(_translate("PAGEANT", "./bin/key"))
        self.i_qr_code_key.setPlaceholderText(_translate("PAGEANT", "Select MAF reference data"))
        self.Function.setTabText(self.Function.indexOf(self.API), _translate("PAGEANT", "API"))


class MyMainForm(QMainWindow, Ui_PAGEANT):
    def __init__(self, parent=None):
        super(MyMainForm, self).__init__(parent)
        self.setupUi(self)
        self.s_vcf.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_vcf, 'input_file', "Select the input file", os.getcwd(),
                       filter='All file formats supported (*.vcf *.vcf.gz *.bed *.txt);;'
                              ' vcf (*.vcf *.vcf.gz);; 23andme (*.txt);; PLINK 1 binary (*.bed)'))
        self.s_config.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_config, 'config_file',
                       "Select the configuration file", os.getcwd()))
        self.s_output.clicked.connect(
            self.open_(QFileDialog.getExistingDirectory, self.i_output, 'output',
                       "Select the output directory", os.getcwd()))
        self.s_database_qual.clicked.connect(main.partial(self.o_database, 'qual'))
        self.s_database_quan.clicked.connect(main.partial(self.o_database, 'quan'))
        self.s_ref_qual.clicked.connect(main.partial(self.o_ref, 'qual'))
        self.s_ref_quan.clicked.connect(main.partial(self.o_ref, 'quan'))
        self.s_qr_dir.clicked.connect(
            self.open_(QFileDialog.getExistingDirectory, self.i_qr_dir, 'qr_dir',
                       "Specify the directory where needed files can be found", os.getcwd()))
        self.s_qr_user.clicked.connect(
            self.open_(QFileDialog.getSaveFileName, self.i_qr_user, 'qr_user',
                       "Specify the path for saving User QR code",
                       os.path.join(os.getcwd(), 'user_qr_code.png'), 'image(*.png)'))
        self.s_qr_snps.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_qr_snps, 'qr_snps',
                       "Select the text which include the needed SNPs list", os.getcwd()))
        self.s_qr_give.clicked.connect(
            self.open_(QFileDialog.getSaveFileName, self.i_qr_give, 'qr_give',
                       "Specify the path for saving SNP list QR code",
                       os.path.join(os.getcwd(), 'snp_qr_code.png'), 'image(*.png)'))
        self.s_maf_sample.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_maf_sample, 'maf_ref',
                       "Select the MAF reference file", os.getcwd(),
                       filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; PLINK 2 binary (*.pgen)'))
        self.s_ps_ref.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_ps_ref, 'ps_ref',
                       "Select the reference file for population stratification analysis", os.getcwd(),
                       filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; PLINK 2 binary (*.pgen)'))
        self.s_ps_pop.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_ps_pop, 'population_file',
                       "Select the population file for population stratification analysis", os.getcwd()))
        self.s_concord_ref.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_concord_ref, 'concord_ref',
                       "Select the concordance analysis reference file", os.getcwd(),
                       filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; PLINK 2 binary (*.pgen)'))
        self.s_otherdb.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_otherdb, 'query_db',
                       "Select the third-party database file", os.getcwd()))
        self.s_umap_ref.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_umap_ref, None,
                       "Select the reference file for PCA & UMAP plot", os.getcwd(),
                       filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; PLINK 2 binary (*.pgen)'))
        self.s_umap_sample.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_umap_sample, None, "Select the sample file", os.getcwd(),
                       filter='All file formats supported (*.vcf *.vcf.gz *.bed *.txt);;'
                              ' vcf (*.vcf *.vcf.gz);; 23andme (*.txt);; PLINK 1 binary (*.bed)'))
        self.s_umap_metadata.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_umap_metadata, None,
                       "Select the metadata for reference genotype data", os.getcwd(),
                       filter='All file formats supported (*.vcf *.vcf.gz *.bed *.txt);;'
                              ' vcf (*.vcf *.vcf.gz);; 23andme (*.txt);; PLINK 1 binary (*.bed)'))
        self.s_add_rsid_file.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_add_rsid_file, None,
                       "Select the GWAS file", os.getcwd()))
        self.s_add_rsid_db.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_add_rsid_db, None,
                       "Select the dbSNP file", os.getcwd()))
        self.s_qr_code_snp.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_qr_code_snp, None,
                       "Select the text which include the needed SNPs list", os.getcwd()))
        self.s_qr_code_key.clicked.connect(
            self.open_(QFileDialog.getOpenFileName, self.i_qr_code_key, None,
                       "Select the file which contains public key", os.getcwd())
        )

        self.b_analyze.clicked.connect(self.run)
        self.b_default_qual.clicked.connect(self.set_default_qual)
        self.b_default_quan.clicked.connect(self.set_default_quan)
        self.b_default_sample_qc.clicked.connect(self.set_default_sample_qc)
        self.b_default_query_db.clicked.connect(self.set_default_query_db)
        self.b_default_qr_code.clicked.connect(self.set_default_qr_code)
        # todo
        self.b_run_umap.clicked.connect(self.run_umap)
        self.b_run_add_rsid.clicked.connect(self.run_add_rsid)
        self.b_run_qr_code.clicked.connect(self.run_qr_code)

        self.timer = QtCore.QBasicTimer()

        self.inter = 50
        self.progress_value = ProgressValue()
        self.progress_now = ProgressNow(self.inter)

        # progress config dict
        self.parameters = default_config

    def o_database(self, otype: str = 'qual' or 'quan'):
        get_dir_path = QFileDialog.getExistingDirectory(self, "Select the database directory", os.getcwd())
        if otype == 'quan':
            self.i_database_quan.setText(get_dir_path)
        elif otype == 'qual':
            self.i_database_qual.setText(get_dir_path)
        self.parameters[f'{otype}_data'] = get_dir_path

    def o_ref(self, otype: str = 'qual' or 'quan'):
        get_dir_path = QFileDialog.getExistingDirectory(self, "Select the reference population directory", os.getcwd())
        if otype == 'quan':
            self.i_ref_quan.setText(get_dir_path)
        elif otype == 'qual':
            self.i_ref_qual.setText(get_dir_path)
        self.parameters[f'{otype}_ref'] = get_dir_path

    def open_raw(self, open_func: main.Callable, dest_obj: QtWidgets.QLineEdit, dest: main.Optional[str] = None,
                 title: main.Optional[str] = 'Open', default: str = '', *args, filter_: str = '', **kwargs) -> None:
        if open_func != QFileDialog.getExistingDirectory:
            get_directory_path = open_func(self, title, default, filter_)
        else:
            get_directory_path = open_func(self, title, default)
        get_path = get_directory_path if type(get_directory_path) == str else get_directory_path[0]
        dest_obj.setText(get_path)
        if dest:
            self.parameters[dest] = get_path

    def open_(self, *args, **kargs) -> main.Callable:
        return main.partial(self.open_raw, *args, **kargs)

    def set_default_qual(self):
        self.i_database_qual.setText("./algorithm/Qualitative")
        self.i_ref_qual.setText("./population_genome")

    def set_default_quan(self):
        self.i_database_quan.setText("./algorithm/Quantitative")
        self.i_ref_quan.setText("./population_genome")
        self.s_quan_pop.setCurrentIndex(1)

    def set_default_sample_qc(self):
        self.i_maf_sample.setText(default_config['maf_ref'])
        self.i_ps_ref.setText(str(default_config['ps_ref']))
        self.i_ps_pop.setText(default_config['population_file'])
        self.i_concord_ref.setText(str(default_config['concord_ref']))
        self.i_pop.setText(default_config['population_col'])
        self.i_pop_id.setText(str(default_config['population_id']))
        self.s_pop_sep.setCurrentIndex(0)

    def set_default_query_db(self):
        self.c_query_database.setChecked(True)
        self.c_clinvar.setChecked(True)
        self.c_pharmgkb.setChecked(True)
        self.c_otherdb.setChecked(True)
        self.i_otherdb.setText('./algorithm/Query_db/Phewas_Catalog_part.tsv')

    def set_default_qr_code(self):
        self.c_qr_code.setChecked(True)
        self.i_qr_snps.setText('./personal_genome/fingerprint_snps.txt')
        self.i_qr_give.setText('./output/qr_code/snp_qr_code.png')
        self.i_qr_dir.setText("./output/qr_code")
        self.i_qr_user.setText("./output/qr_code/user_qr_code.png")

    def timerEvent(self, a0: 'QtCore.QTimerEvent') -> None:
        if self.progress_now.value() < self.progress_value.value:
            self.progress_now.set(self.progress_value.value)
        elif self.progress_now.value() < self.progress_value.value + 4:
            self.progress_now.add()
        if self.progress.value() != self.progress_now.value():
            self.progress.setValue(self.progress_now.value())

    def cancel_failed(self):
        QtWidgets.QMessageBox.critical(self, 'Error', 'Cancel function disable!', QtWidgets.QMessageBox.Close)

    def report(self):
        self.b_analyze.setEnabled(True)
        self.progress.setWindowFlags(QtCore.Qt.WindowMinimizeButtonHint |
                                     QtCore.Qt.WindowCloseButtonHint)
        self.thread.deleteLater()
        self.timer.stop()
        if self.thread.statue:
            QtWidgets.QMessageBox.information(self, 'Result', self.thread.res)
        else:
            QtWidgets.QMessageBox.critical(self, 'Error', self.thread.res, QtWidgets.QMessageBox.Close)
        reply = QtWidgets.QMessageBox.question(self, 'Question', 'Do you want to open report now?',
                                               QtWidgets.QMessageBox.Ok |
                                               QtWidgets.QMessageBox.Close,
                                               QtWidgets.QMessageBox.Ok)
        self.b_analyze.setEnabled(True)
        self.thread.deleteLater()
        if reply == QtWidgets.QMessageBox.Ok:
            if main.platform == 'Windows' or main.platform == 'Linux':
                webbrowser.open_new(os.path.join(self.parameters['output'], 'genetic_report', 'Report.html'))
            elif main.platform == 'Darwin':
                webbrowser.get().open('file:///' + os.path.abspath(os.path.join(self.parameters['output'],
                                                                                'genetic_report', 'Report.html')))
            else:
                QtWidgets.QMessageBox.warning(self, 'Failed', 'Unsupported operating system to open file',
                                              QtWidgets.QMessageBox.Ok)
        else:
            return

    def run(self):
        self.parameters.update(dict(zip(['name', 'input_file', 'config_file', 'output', 'qual_data', 'qual_ref',
                                         'quan_data', 'quan_ref', 'sample_qc', 'vep', 'maf_ref', 'ps_ref',
                                         'concord_ref', 'population_col', 'population_file',
                                         'population_id', 'query_database', 'pharmgkb', 'clinvar', 'query_db',
                                         'qr_code', 'qr_dir', 'qr_snps', 'qr_give', 'qr_user', 'quan_pop'],
                                        [self.i_name.text(), self.i_vcf.text(), self.i_config.text(),
                                         self.i_output.text(), self.i_database_qual.text(), self.i_ref_qual.text(),
                                         self.i_database_quan.text(), self.i_ref_quan.text(),
                                         self.c_sample_qc.isChecked(), self.c_vep.isChecked(),
                                         self.i_maf_sample.text(), self.i_ps_ref.text(), self.i_concord_ref.text(),
                                         self.i_pop.text(), self.i_ps_pop.text(), self.i_pop_id.text(),
                                         self.c_query_database.isChecked(), self.c_pharmgkb.isChecked(),
                                         self.c_clinvar.isChecked(),
                                         self.i_otherdb.text() if self.c_otherdb.isChecked() else '',
                                         self.c_qr_code.isChecked(), self.i_qr_dir.text(), self.i_qr_snps.text(),
                                         self.i_qr_give.text(), self.i_qr_user.text(), self.s_quan_pop.currentText()
                                         ])))
        self.progress_value = ProgressValue()
        self.progress_now = ProgressNow(self.inter)
        self.progress = QtWidgets.QProgressDialog('Analysis start.', 'Cancel', 0, 100)
        self.progress.setWindowModality(QtCore.Qt.NonModal)
        self.progress.setAutoClose(True)
        self.progress.setAutoReset(False)
        self.progress.setWindowFlags(QtCore.Qt.WindowMinimizeButtonHint)
        # todo: How to cancel the progress?
        self.progress.canceled.connect(self.cancel_failed)
        self.progress.setFixedSize(500, 150)
        self.progress.show()
        self.b_analyze.setEnabled(False)

        handler = GUIHandler(self.progress, self.progress_value)
        self.thread = MyThread(handler, self.parameters)

        self.timer.start(self.inter, self)
        self.thread.start()
        self.thread.finished.connect(self.report)

    def run_umap(self):
        temp_dir = main.mkdtemp(suffix='umap')
        self.thread = API_Thread(umap, self.i_umap_ref.text(), self.i_umap_sample.text(),
                                 self.i_umap_metadata.text(), temp_dir, temp_dir, os.getcwd(),
                                 prune=True, pop_col='population', pop_id='IID')
        self.setEnabled(False)
        self.thread.start()
        self.thread.finished.connect(main.partial(self.api_finish, temp_dir))

    def api_finish(self, temp_dir: main.Optional[str] = None):
        self.setEnabled(True)
        QtWidgets.QMessageBox.information(self, 'Finish', self.thread.res)
        if temp_dir:
            main.rm_dir(temp_dir)

    def run_add_rsid(self):
        import src.add_rsid as add_rsid
        self.thread = API_Thread(add_rsid.run, self.i_add_rsid_file.text(), self.i_add_rsid_db.text(),
                                 os.path.join(os.getcwd(), 'annotated.tsv.gz'),
                                 end_='Finish! Result (annotated.tsv.gz) has been saved in the working directory.')
        self.setEnabled(False)
        self.thread.start()
        self.thread.finished.connect(self.api_finish)

    def run_qr_code(self):
        import src.qr_code as crypto
        self.thread = API_Thread(crypto.request, self.i_qr_code_key.text(), self.i_qr_code_snp.text(),
                                 os.path.join(os.getcwd(), 'SNP_QR_CODE.png'), None)
        self.setEnabled(False)
        self.thread.start()
        self.thread.finished.connect(self.api_finish)


if __name__ == "__main__":
    freeze_support()
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = MyMainForm()
    MainWindow.show()
    sys.exit(app.exec_())
