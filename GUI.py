# -*- coding: utf-8 -*-
import os
import sys
import webbrowser
import logging
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5 import __file__ as qt_file
from PyQt5.QtWidgets import QMainWindow, QFileDialog
from multiprocessing import freeze_support
from time import sleep
import pageant as main


main.plt.switch_backend('Agg')
if main.platform == 'Darwin':
    os.chdir(main.raw_dir)
y_offset = 0 if main.platform == 'Windows' else 10
QtGui.QGuiApplication.addLibraryPath(os.path.join(os.path.dirname(qt_file), 'Qt', 'plugins', 'platforms'))
QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)

default_config = {value: key[value] for key in main.load_config().values() for value in key}
default_config.update({
    'name': 'Test',
    'input_file': './personal_genome/sample.vcf.gz',
    'config_file': './bin/config.ini',
    'quan_data': './algorithm_database/Quantitative',
    'qual_data': './algorithm_database/Qualitative',
    'quan_ref': './population_genome',
    'qual_ref': './population_genome',
    # 'qr_snps_txt': './personal_genome/fingerprint_snps.txt',
    'maf_ref': './personal_genome/hapmap3.vcf.gz',
    'ps_ref': './personal_genome/hapmap3.vcf.gz',
    'concord_ref': './personal_genome/concordance.vcf.gz'
})


def index_to_sep(ind: int) -> str:
    return {0: '\t', 1: ',', 2: ' '}[ind]


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
        font.setPointSize(11)
        font.setBold(True)
        font1 = QtGui.QFont()
        font1.setFamily("Calibri")
        font1.setPointSize(10)
        font1.setBold(False)
        font2 = QtGui.QFont()
        font2.setFamily("Calibri")
        font2.setPointSize(11)
        font2.setBold(False)
        font2.setKerning(True)
        font2.setStyleStrategy(QtGui.QFont.PreferDefault)
        font3 = QtGui.QFont()
        font3.setFamily("Calibri")
        font3.setPointSize(11)
        font3.setBold(False)
        font4 = QtGui.QFont()
        font4.setFamily("Calibri")
        font4.setPointSize(11)
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
        self.l_database_quan.setGeometry(QtCore.QRect(30, 80, 201, 25))

        self.l_database_quan.setFont(font)
        self.l_database_quan.setObjectName("l_database_quan")
        self.s_ref_quan = QtWidgets.QToolButton(self.g_quan)
        self.s_ref_quan.setGeometry(QtCore.QRect(430, 230, 41, 27))

        self.s_ref_quan.setFont(font1)
        self.s_ref_quan.setObjectName("s_ref_quan")
        self.i_database_quan = QtWidgets.QLineEdit(self.g_quan)
        self.i_database_quan.setGeometry(QtCore.QRect(69, 111, 365, 25))

        self.i_database_quan.setFont(font1)
        self.i_database_quan.setObjectName("i_database_quan")
        self.i_ref_quan = QtWidgets.QLineEdit(self.g_quan)
        self.i_ref_quan.setGeometry(QtCore.QRect(69, 231, 365, 25))

        self.i_ref_quan.setFont(font1)
        self.i_ref_quan.setObjectName("i_ref_quan")
        self.s_database_quan = QtWidgets.QToolButton(self.g_quan)
        self.s_database_quan.setGeometry(QtCore.QRect(430, 110, 41, 27))

        self.s_database_quan.setFont(font1)
        self.s_database_quan.setObjectName("s_database_quan")
        self.l_ref_quan = QtWidgets.QLabel(self.g_quan)
        self.l_ref_quan.setGeometry(QtCore.QRect(30, 200, 201, 25))

        self.l_ref_quan.setFont(font)
        self.l_ref_quan.setObjectName("l_ref_quan")
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

        self.l_qr_user = QtWidgets.QLabel(self.c_qr_code)
        self.l_qr_user.setGeometry(QtCore.QRect(30, 210, 281, 25))

        self.l_qr_user.setFont(font)
        self.l_qr_user.setObjectName("l_qr_user")
        self.i_qr_user = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_qr_user.setGeometry(QtCore.QRect(50, 241, 365, 25))

        self.i_qr_user.setFont(font1)
        self.i_qr_user.setObjectName("i_qr_user")
        self.s_qr_user = QtWidgets.QToolButton(self.c_qr_code)
        self.s_qr_user.setGeometry(QtCore.QRect(410, 240, 41, 27))

        self.s_qr_user.setFont(font1)
        self.s_qr_user.setObjectName("s_qr_user")
        self.l_qr_dr = QtWidgets.QLabel(self.c_qr_code)
        self.l_qr_dr.setGeometry(QtCore.QRect(30, 150, 281, 25))

        self.l_qr_dr.setFont(font)
        self.l_qr_dr.setObjectName("l_qr_dr")
        self.i_qr_dr = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_qr_dr.setGeometry(QtCore.QRect(50, 181, 365, 25))

        self.i_qr_dr.setFont(font1)
        self.i_qr_dr.setObjectName("i_qr_dr")
        self.s_qr_dr = QtWidgets.QToolButton(self.c_qr_code)
        self.s_qr_dr.setGeometry(QtCore.QRect(410, 180, 41, 27))

        self.s_qr_dr.setFont(font1)
        self.s_qr_dr.setObjectName("s_qr_dr")
        self.i_key = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_key.setGeometry(QtCore.QRect(50, 61, 361, 25))

        self.i_key.setFont(font1)
        self.i_key.setObjectName("i_key")
        self.s_key = QtWidgets.QToolButton(self.c_qr_code)
        self.s_key.setGeometry(QtCore.QRect(407, 60, 41, 27))

        self.s_key.setFont(font1)
        self.s_key.setObjectName("s_key")
        self.l_key = QtWidgets.QLabel(self.c_qr_code)
        self.l_key.setGeometry(QtCore.QRect(30, 30, 291, 25))

        self.l_key.setFont(font)
        self.l_key.setObjectName("l_key")
        self.l_qr_dir = QtWidgets.QLabel(self.c_qr_code)
        self.l_qr_dir.setGeometry(QtCore.QRect(30, 270, 251, 25))

        self.l_qr_dir.setFont(font)
        self.l_qr_dir.setObjectName("l_qr_dir")
        self.i_qr_dir = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_qr_dir.setGeometry(QtCore.QRect(53, 301, 361, 25))

        self.i_qr_dir.setFont(font1)
        self.i_qr_dir.setObjectName("i_qr_dir")
        self.s_qr_dir = QtWidgets.QToolButton(self.c_qr_code)
        self.s_qr_dir.setGeometry(QtCore.QRect(410, 300, 41, 27))

        self.s_qr_dir.setFont(font1)
        self.s_qr_dir.setObjectName("s_qr_dir")
        self.l_qr_snps = QtWidgets.QLabel(self.c_qr_code)
        self.l_qr_snps.setGeometry(QtCore.QRect(30, 90, 281, 25))

        self.l_qr_snps.setFont(font)
        self.l_qr_snps.setObjectName("l_qr_snps")
        self.i_qr_snps = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_qr_snps.setGeometry(QtCore.QRect(50, 121, 361, 25))

        self.i_qr_snps.setFont(font1)
        self.i_qr_snps.setObjectName("i_qr_snps")
        self.s_qr_snps = QtWidgets.QToolButton(self.c_qr_code)
        self.s_qr_snps.setGeometry(QtCore.QRect(406, 120, 41, 27))

        self.s_qr_snps.setFont(font1)
        self.s_qr_snps.setObjectName("s_qr_snps")

        self.b_default_qr_code = QtWidgets.QPushButton(self.QR_code)
        self.b_default_qr_code.setGeometry(QtCore.QRect(230, 380, 91, 31))

        self.b_default_qr_code.setFont(font)
        self.b_default_qr_code.setObjectName("b_default_qr_code")
        self.Function.addTab(self.QR_code, "")
        PAGEANT.setCentralWidget(self.MainWindow)

        self.retranslateUi(PAGEANT)
        self.Function.setCurrentIndex(0)
        self.s_pop_sep.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(PAGEANT)


    def retranslateUi(self, PAGEANT):
        _translate = QtCore.QCoreApplication.translate
        PAGEANT.setWindowTitle(_translate("PAGEANT", f"PAGEANT ({main.version})"))
        self.i_vcf.setToolTip(_translate("PAGEANT", "Select the genotype file that you want to analyze"))
        self.i_vcf.setText(_translate("PAGEANT", "./personal_genome/sample.vcf.gz"))
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
        self.i_name.setText(_translate("PAGEANT", "test"))
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
        self.i_database_qual.setText(_translate("PAGEANT", "./algorithm_database/Qualitative"))
        self.i_database_qual.setPlaceholderText(_translate("PAGEANT", "Specify database for qualitative traits"))
        self.i_ref_qual.setToolTip(_translate("PAGEANT", "Specify directory of reference for qualitative traits"))
        self.i_ref_qual.setText(_translate("PAGEANT", "./population_genome"))
        self.i_ref_qual.setPlaceholderText(_translate("PAGEANT", "Specify reference for qualitative traits"))
        self.s_database_qual.setText(_translate("PAGEANT", "..."))
        self.l_ref_qual.setText(_translate("PAGEANT", "Reference Population directory"))
        self.Function.setTabText(self.Function.indexOf(self.Qualitative), _translate("PAGEANT", "Qualitative"))
        self.b_default_quan.setToolTip(_translate("PAGEANT", "Set default values"))
        self.b_default_quan.setText(_translate("PAGEANT", "Default"))
        self.g_quan.setTitle(_translate("PAGEANT", "Quantitative traits"))
        self.l_database_quan.setText(_translate("PAGEANT", "Database directory"))
        self.s_ref_quan.setText(_translate("PAGEANT", "..."))
        self.i_database_quan.setToolTip(_translate("PAGEANT", "Specify directory of database for quantitative traits"))
        self.i_database_quan.setText(_translate("PAGEANT", "./algorithm_database/Quantitative"))
        self.i_database_quan.setPlaceholderText(_translate("PAGEANT", "Specify database for quantitative traits"))
        self.i_ref_quan.setToolTip(_translate("PAGEANT", "Specify directory of reference for quantitative traits"))
        self.i_ref_quan.setText(_translate("PAGEANT", "./population_genome"))
        self.i_ref_quan.setPlaceholderText(_translate("PAGEANT", "Specify reference for quantitative traits"))
        self.s_database_quan.setText(_translate("PAGEANT", "..."))
        self.l_ref_quan.setText(_translate("PAGEANT", "Reference Population directory"))
        self.Function.setTabText(self.Function.indexOf(self.Quantitative), _translate("PAGEANT", "Quantitative"))

        # self.c_sample_qc.setToolTip(_translate("PAGEANT", "Activate function of sample QC"))
        self.c_sample_qc.setTitle(_translate("PAGEANT", "Sample QC"))
        self.c_vep.setToolTip(_translate("PAGEANT", "Using VEP to inspect the variants\' type in sample data"))
        self.c_vep.setText(_translate("PAGEANT", "Use VEP (Slow with too many variants or bad network!)"))
        self.PS.setTitle(_translate("PAGEANT", "Population stratification analysis"))
        self.l_ps_ref.setText(_translate("PAGEANT", "Reference data"))
        self.i_ps_ref.setToolTip(_translate("PAGEANT", "Specify the reference data for ps"))
        self.i_ps_ref.setText(_translate("PAGEANT", "./personal_genome/hapmap3.vcf.gz"))
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
        self.i_ps_pop.setText(_translate("PAGEANT", "./personal_genome/hapmap3_samples.txt"))
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
        self.i_maf_sample.setText(_translate("PAGEANT", "./personal_genome/hapmap3.vcf.gz"))
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
        self.i_otherdb.setText(_translate("PAGEANT", "./algorithm_database/Query_database/Phewas_catalog_part.tsv"))
        self.i_otherdb.setPlaceholderText(_translate("PAGEANT", "Select third-party database file"))
        self.s_otherdb.setText(_translate("PAGEANT", "..."))
        self.Function.setTabText(self.Function.indexOf(self.Query_database), _translate("PAGEANT", "Query DB"))

        # self.c_qr_code.setToolTip(_translate("PAGEANT", "Activate function of generating QR code"))
        self.c_qr_code.setTitle(_translate("PAGEANT", "Generate QR code"))
        self.l_qr_user.setText(_translate("PAGEANT", "Directory for user QR code file"))
        self.i_qr_user.setToolTip(_translate("PAGEANT", "Specify the text which include the needed SNPs list"))
        self.i_qr_user.setText(_translate("PAGEANT", "./output/qr_code/user"))
        self.s_qr_user.setText(_translate("PAGEANT", "..."))
        self.l_qr_dr.setText(_translate("PAGEANT", "Directory for doctor QR code file"))
        self.i_qr_dr.setToolTip(_translate("PAGEANT", "Specify the text which include the needed SNPs list"))
        self.i_qr_dr.setText(_translate("PAGEANT", "./output/qr_code/doctor"))
        self.s_qr_dr.setText(_translate("PAGEANT", "..."))
        self.i_key.setToolTip(_translate("PAGEANT", "Specify your key file"))
        self.i_key.setText(_translate("PAGEANT", "./output/qr_code/doctor/key"))
        self.s_key.setText(_translate("PAGEANT", "..."))
        self.l_key.setText(_translate("PAGEANT", "Doctor\'s public/private Key file"))
        self.l_qr_dir.setText(_translate("PAGEANT", "Directory for decrypted user genotype"))
        self.i_qr_dir.setToolTip(_translate("PAGEANT", "Specify directory where the result save"))
        self.i_qr_dir.setText(_translate("PAGEANT", "./output/qr_code/doctor"))
        self.s_qr_dir.setText(_translate("PAGEANT", "..."))
        self.l_qr_snps.setText(_translate("PAGEANT", "Doctor\'s SNP list file"))
        self.i_qr_snps.setToolTip(_translate("PAGEANT", "Specify the text which include the needed SNPs list"))
        self.i_qr_snps.setText(_translate("PAGEANT", "./personal_genome/fingerprint_snps.txt"))
        self.s_qr_snps.setText(_translate("PAGEANT", "..."))

        self.b_default_qr_code.setToolTip(_translate("PAGEANT", "Set default values"))
        self.b_default_qr_code.setText(_translate("PAGEANT", "Default"))
        self.Function.setTabText(self.Function.indexOf(self.QR_code), _translate("PAGEANT", "QR code"))


class MyMainForm(QMainWindow, Ui_PAGEANT):
    def __init__(self, parent=None):
        super(MyMainForm, self).__init__(parent)
        self.setupUi(self)

        self.s_vcf.clicked.connect(self.o_vcf)
        self.s_config.clicked.connect(self.o_config)
        self.s_output.clicked.connect(self.o_output)
        self.s_database_qual.clicked.connect(main.partial(self.o_database, 'qual'))
        self.s_database_quan.clicked.connect(main.partial(self.o_database, 'quan'))
        self.s_ref_qual.clicked.connect(main.partial(self.o_ref, 'qual'))
        self.s_ref_quan.clicked.connect(main.partial(self.o_ref, 'quan'))
        self.s_qr_dr.clicked.connect(self.o_qr_dr)
        self.s_qr_user.clicked.connect(self.o_qr_user)
        self.s_qr_dir.clicked.connect(self.o_qr_dir)
        self.s_qr_snps.clicked.connect(self.o_qr_snps)
        self.s_key.clicked.connect(self.o_qr_key)
        self.s_maf_sample.clicked.connect(self.o_maf)
        self.s_ps_ref.clicked.connect(self.o_ps)
        self.s_ps_pop.clicked.connect(self.o_text)
        self.s_concord_ref.clicked.connect(self.o_concord)
        self.s_otherdb.clicked.connect(self.o_otherdb)
        self.b_analyze.clicked.connect(self.run)
        self.b_default_qual.clicked.connect(self.set_default_qual)
        self.b_default_quan.clicked.connect(self.set_default_quan)
        self.b_default_sample_qc.clicked.connect(self.set_default_sample_qc)
        self.b_default_query_db.clicked.connect(self.set_default_query_db)
        self.b_default_qr_code.clicked.connect(self.set_default_qr_code)

        self.timer = QtCore.QBasicTimer()

        self.inter = 50
        self.progress_value = ProgressValue()
        self.progress_now = ProgressNow(self.inter)

        # progress config dict
        self.parameters = default_config

    def o_vcf(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the input file", os.getcwd(),
                                                         filter='All file formats supported '
                                                                '(*.vcf *.vcf.gz *.bed *.txt);; '
                                                                'vcf (*.vcf *.vcf.gz);; 23andme (*.txt);; '
                                                                'PLINK 1 binary (*.bed)')
        self.i_vcf.setText(get_directory_path[0])
        self.parameters['input_file'] = get_directory_path[0]

    def o_config(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the configuration file", os.getcwd())
        self.i_config.setText(get_directory_path[0])
        self.parameters['config_file'] = get_directory_path[0]

    def o_output(self):
        get_dir_path = QFileDialog.getExistingDirectory(self, "Select the output directory", os.getcwd())
        self.i_output.setText(get_dir_path)
        self.parameters['output'] = get_dir_path

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

    def o_maf(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the MAF reference file", os.getcwd(),
                                                         filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; '
                                                                'PLINK 2 binary (*.pgen)')
        self.i_maf_sample.setText(get_directory_path[0])
        self.parameters['maf_ref'] = get_directory_path[0]

    def o_ps(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the reference file for "
                                                               "population stratification analysis", os.getcwd(),
                                                         filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; '
                                                                'PLINK 2 binary (*.pgen)')
        self.i_ps_ref.setText(get_directory_path[0])
        self.parameters['ps_ref'] = get_directory_path[0]

    def o_concord(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the concordance analysis reference file", os.getcwd(),
                                                         filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; '
                                                                'PLINK 2 binary (*.pgen)')
        self.i_concord_ref.setText(get_directory_path[0])
        self.parameters['concord_ref'] = get_directory_path[0]

    def o_qr_key(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the doctor's key file", os.getcwd())
        self.i_key.setText(get_directory_path[0])

    def o_qr_dir(self):
        get_dir_path = QFileDialog.getExistingDirectory(self, "Select the directory where the result saves",
                                                        os.getcwd())
        self.i_qr_dir.setText(get_dir_path)

    def o_qr_snps(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the text which include the needed SNPs list",
                                                         os.getcwd())
        self.i_qr_snps.setText(get_directory_path[0])

    def o_qr_dr(self):
        get_dir_path = QFileDialog.getExistingDirectory(self, "Select the directory for doctor QR code",
                                                        os.getcwd())
        self.i_qr_dr.setText(get_dir_path)
        self.parameters['qr_dr'] = get_dir_path

    def o_qr_user(self):
        get_dir_path = QFileDialog.getExistingDirectory(self, "Select the directory for user QR code",
                                                        os.getcwd())
        self.i_qr_user.setText(get_dir_path)
        self.parameters['qr_user'] = get_dir_path

    def o_otherdb(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the third-party database file", os.getcwd())
        self.i_otherdb.setText(get_directory_path[0])
        self.parameters['query_db'] = get_directory_path[0]

    def o_text(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the population file for "
                                                               "population stratification analysis", os.getcwd())
        self.i_ps_pop.setText(get_directory_path[0])
        self.parameters['population_file'] = get_directory_path[0]

    def set_default_qual(self):
        self.i_database_qual.setText("./algorithm_database/Qualitative")
        self.i_ref_qual.setText("./population_genome")

    def set_default_quan(self):
        self.i_database_quan.setText("./algorithm_database/Quantitative")
        self.i_ref_quan.setText("./population_genome")

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
        self.i_otherdb.setText('./algorithm_database/Query_database/Phewas_catalog_part.tsv')

    def set_default_qr_code(self):
        self.c_qr_code.setChecked(True)
        self.i_key.setText('./output/qr_code/doctor/key')
        self.i_qr_snps.setText('./personal_genome/fingerprint_snps.txt')
        self.i_qr_dir.setText('./output/qr_code/doctor')
        self.i_qr_dr.setText("./output/qr_code/doctor")
        self.i_qr_user.setText("./output/qr_code/user")

    def timerEvent(self, a0: 'QtCore.QTimerEvent') -> None:
        if self.progress_now.value() < self.progress_value.value:
            self.progress_now.set(self.progress_value.value)
        elif self.progress_now.value() < self.progress_value.value + 4:
            self.progress_now.add()
        if self.progress.value() != self.progress_now.value():
            self.progress.setValue(self.progress_now.value())

    def cancel_failed(self):
        QtWidgets.QMessageBox.critical(self, 'Error', 'Cancel function disable!',
                                       QtWidgets.QMessageBox.Close)

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
                                         'quan_data',  'quan_ref', 'sample_qc', 'vep', 'maf_ref', 'ps_ref',
                                         'concord_ref', 'population_col', 'population_file',
                                         'population_id', 'query_database', 'pharmgkb', 'clinvar', 'query_db',
                                         'qr_code', 'qr_key', 'qr_snps', 'qr_dr', 'qr_user', 'qr_dir'],
                                        [self.i_name.text(), self.i_vcf.text(), self.i_config.text(),
                                         self.i_output.text(), self.i_database_qual.text(), self.i_ref_qual.text(),
                                         self.i_database_quan.text(), self.i_ref_quan.text(),
                                         self.c_sample_qc.isChecked(),  self.c_vep.isChecked(),
                                         self.i_maf_sample.text(), self.i_ps_ref.text(), self.i_concord_ref.text(),
                                         self.i_pop.text(), self.i_ps_pop.text(), self.i_pop_id.text(),
                                         self.c_query_database.isChecked(), self.c_pharmgkb.isChecked(),
                                         self.c_clinvar.isChecked(),
                                         self.i_otherdb.text() if self.c_otherdb.isChecked() else '',
                                         self.c_qr_code.isChecked(), self.i_key.text(), self.i_qr_snps.text(),
                                         self.i_qr_dr.text(), self.i_qr_user.text(), self.i_qr_dir.text()
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


if __name__ == "__main__":
    freeze_support()
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = MyMainForm()
    MainWindow.show()
    sys.exit(app.exec_())