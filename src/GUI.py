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
import main


main.plt.switch_backend('Agg')
if main.platform == 'Darwin':
    os.chdir(main.raw_dir)
QtGui.QGuiApplication.addLibraryPath(os.path.join(os.path.dirname(qt_file), 'Qt', 'plugins', 'platforms'))
QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)

default_config = main.load_config()['parameter']
default_config.update({
    'name': 'Test',
    'input_file': './personal_genome/sample.vcf.gz',
    'data_dir': './algorithm_database',
    'ref': './population_genome',
    'qr_snps_txt': './personal_genome/fingerprint_snps.txt',
    'maf_ref': './personal_genome/hapmap3.vcf.gz',
    'pca_ref': './personal_genome/hapmap3.vcf.gz',
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
        font3 = QtGui.QFont()
        font3.setFamily("Calibri")
        font3.setPointSize(10)
        font3.setBold(True)
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
        self.Basic = QtWidgets.QWidget()
        self.Basic.setObjectName(u"Basic")
        self.Basic.setFont(font)
        self.i_vcf = QtWidgets.QLineEdit(self.Basic)
        self.i_vcf.setObjectName("i_vcf")
        self.i_vcf.setGeometry(QtCore.QRect(79, 111, 365, 25))
        self.i_vcf.setFont(font1)
        self.b_analyze = QtWidgets.QPushButton(self.Basic)
        self.b_analyze.setObjectName("b_analyze")
        self.b_analyze.setGeometry(QtCore.QRect(230, 380, 91, 31))
        self.b_analyze.setFont(font)
        self.s_vcf = QtWidgets.QToolButton(self.Basic)
        self.s_vcf.setObjectName("s_vcf")
        self.s_vcf.setGeometry(QtCore.QRect(440, 110, 41, 27))
        self.s_vcf.setFont(font1)
        self.l_vcf = QtWidgets.QLabel(self.Basic)
        self.l_vcf.setObjectName("l_vcf")
        self.l_vcf.setGeometry(QtCore.QRect(40, 80, 121, 20))
        self.l_vcf.setFont(font)
        self.l_output = QtWidgets.QLabel(self.Basic)
        self.l_output.setObjectName("l_output")
        self.l_output.setGeometry(QtCore.QRect(40, 290, 171, 20))
        self.l_output.setFont(font)
        self.i_output = QtWidgets.QLineEdit(self.Basic)
        self.i_output.setObjectName("i_output")
        self.i_output.setGeometry(QtCore.QRect(79, 321, 401, 25))
        self.i_output.setFont(font1)
        self.s_output = QtWidgets.QToolButton(self.Basic)
        self.s_output.setObjectName("s_ouput")
        self.s_output.setGeometry(QtCore.QRect(440, 320, 41, 27))
        self.s_output.setFont(font1)
        self.l_name = QtWidgets.QLabel(self.Basic)
        self.l_name.setObjectName("l_name")
        self.l_name.setGeometry(QtCore.QRect(40, 30, 61, 20))
        self.l_name.setFont(font)
        self.i_name = QtWidgets.QLineEdit(self.Basic)
        self.i_name.setObjectName("i_name")
        self.i_name.setGeometry(QtCore.QRect(100, 30, 121, 25))
        self.i_name.setFont(font1)
        self.l_database = QtWidgets.QLabel(self.Basic)
        self.l_database.setObjectName("l_database")
        self.l_database.setGeometry(QtCore.QRect(40, 150, 201, 25))
        self.l_database.setFont(font)
        self.i_database = QtWidgets.QLineEdit(self.Basic)
        self.i_database.setObjectName("i_database")
        self.i_database.setGeometry(QtCore.QRect(79, 181, 365, 25))
        self.i_database.setFont(font1)
        self.s_database = QtWidgets.QToolButton(self.Basic)
        self.s_database.setObjectName("s_database")
        self.s_database.setGeometry(QtCore.QRect(440, 180, 41, 27))
        self.s_database.setFont(font1)
        self.l_ref = QtWidgets.QLabel(self.Basic)
        self.l_ref.setObjectName("l_ref")
        self.l_ref.setGeometry(QtCore.QRect(40, 220, 201, 25))
        self.l_ref.setFont(font)
        self.s_ref = QtWidgets.QToolButton(self.Basic)
        self.s_ref.setObjectName("s_ref")
        self.s_ref.setGeometry(QtCore.QRect(440, 250, 41, 27))
        self.s_ref.setFont(font1)
        self.i_ref = QtWidgets.QLineEdit(self.Basic)
        self.i_ref.setObjectName("i_ref")
        self.i_ref.setGeometry(QtCore.QRect(79, 251, 365, 25))
        self.i_ref.setFont(font1)
        self.Function.addTab(self.Basic, "")
        self.i_output.raise_()
        self.i_vcf.raise_()
        self.b_analyze.raise_()
        self.l_vcf.raise_()
        self.l_output.raise_()
        self.s_output.raise_()
        self.l_name.raise_()
        self.i_name.raise_()
        self.s_vcf.raise_()
        self.l_database.raise_()
        self.i_database.raise_()
        self.s_database.raise_()
        self.l_ref.raise_()
        self.i_ref.raise_()
        self.s_ref.raise_()
        self.QRS = QtWidgets.QWidget()
        self.QRS.setFont(font)
        self.QRS.setObjectName("QRS")
        self.b_default_prs = QtWidgets.QPushButton(self.QRS)
        self.b_default_prs.setGeometry(QtCore.QRect(240, 370, 91, 31))
        self.b_default_prs.setFont(font)
        self.b_default_prs.setObjectName("b_default_prs")
        self.PRS = QtWidgets.QGroupBox(self.QRS)
        self.PRS.setGeometry(QtCore.QRect(30, 140, 501, 211))
        self.PRS.setFont(font)
        self.PRS.setFlat(False)
        self.PRS.setCheckable(False)
        self.PRS.setObjectName("PRS")
        self.l_p_threshold = QtWidgets.QLabel(self.PRS)
        self.l_p_threshold.setGeometry(QtCore.QRect(30, 100, 101, 25))
        self.l_p_threshold.setFont(font2)
        self.l_p_threshold.setObjectName("l_p_threshold")
        self.l_ref_2 = QtWidgets.QLabel(self.PRS)
        self.l_ref_2.setGeometry(QtCore.QRect(30, 140, 281, 25))
        self.l_ref_2.setFont(font2)
        self.l_ref_2.setObjectName("l_ref_2")
        self.l_EA = QtWidgets.QLabel(self.PRS)
        self.l_EA.setGeometry(QtCore.QRect(270, 30, 50, 25))
        self.l_EA.setFont(font2)
        self.l_EA.setObjectName("l_EA")
        self.i_ld_ref = QtWidgets.QLineEdit(self.PRS)
        self.i_ld_ref.setGeometry(QtCore.QRect(51, 170, 351, 25))
        self.i_ld_ref.setFont(font1)
        self.i_ld_ref.setObjectName("i_ld_ref")
        self.l_BETA = QtWidgets.QLabel(self.PRS)
        self.l_BETA.setGeometry(QtCore.QRect(30, 60, 61, 25))
        self.l_BETA.setFont(font2)
        self.l_BETA.setObjectName("l_BETA")
        self.i_p_threshold = QtWidgets.QLineEdit(self.PRS)
        self.i_p_threshold.setGeometry(QtCore.QRect(150, 100, 131, 25))
        self.i_p_threshold.setFont(font1)
        self.i_p_threshold.setObjectName("i_p_threshold")
        self.s_ld_ref_2 = QtWidgets.QToolButton(self.PRS)
        self.s_ld_ref_2.setGeometry(QtCore.QRect(400, 169, 41, 27))

        self.s_ld_ref_2.setFont(font3)
        self.s_ld_ref_2.setObjectName("s_ld_ref_2")
        self.i_EA = QtWidgets.QLineEdit(self.PRS)
        self.i_EA.setGeometry(QtCore.QRect(310, 30, 131, 25))

        self.i_EA.setFont(font1)
        self.i_EA.setObjectName("i_EA")
        self.i_SNP = QtWidgets.QLineEdit(self.PRS)
        self.i_SNP.setGeometry(QtCore.QRect(90, 30, 131, 25))

        self.i_SNP.setFont(font1)
        self.i_SNP.setObjectName("i_SNP")
        self.l_P = QtWidgets.QLabel(self.PRS)
        self.l_P.setGeometry(QtCore.QRect(270, 60, 50, 25))

        self.l_P.setFont(font2)
        self.l_P.setObjectName("l_P")
        self.i_BETA = QtWidgets.QLineEdit(self.PRS)
        self.i_BETA.setGeometry(QtCore.QRect(90, 60, 131, 25))

        self.i_BETA.setFont(font1)
        self.i_BETA.setObjectName("i_BETA")
        self.i_P = QtWidgets.QLineEdit(self.PRS)
        self.i_P.setGeometry(QtCore.QRect(310, 60, 131, 25))

        self.i_P.setFont(font1)
        self.i_P.setObjectName("i_P")
        self.l_SNP = QtWidgets.QLabel(self.PRS)
        self.l_SNP.setGeometry(QtCore.QRect(30, 30, 61, 25))

        self.l_SNP.setFont(font2)
        self.l_SNP.setObjectName("l_SNP")
        self.Clump = QtWidgets.QGroupBox(self.QRS)
        self.Clump.setEnabled(True)
        self.Clump.setGeometry(QtCore.QRect(30, 30, 501, 101))
        self.Clump.setFont(font)
        self.Clump.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.Clump.setObjectName("Clump")
        self.l_p1 = QtWidgets.QLabel(self.Clump)
        self.l_p1.setGeometry(QtCore.QRect(30, 30, 50, 25))

        self.l_p1.setFont(font2)
        self.l_p1.setObjectName("l_p1")
        self.l_p2 = QtWidgets.QLabel(self.Clump)
        self.l_p2.setGeometry(QtCore.QRect(30, 60, 50, 25))

        self.l_p2.setFont(font2)
        self.l_p2.setObjectName("l_p2")
        self.i_p1 = QtWidgets.QLineEdit(self.Clump)
        self.i_p1.setGeometry(QtCore.QRect(80, 30, 131, 25))

        self.i_p1.setFont(font1)
        self.i_p1.setAccessibleDescription("")
        self.i_p1.setObjectName("i_p1")
        self.i_p2 = QtWidgets.QLineEdit(self.Clump)
        self.i_p2.setGeometry(QtCore.QRect(80, 60, 131, 25))

        self.i_p2.setFont(font1)
        self.i_p2.setObjectName("i_p2")
        self.l_r2 = QtWidgets.QLabel(self.Clump)
        self.l_r2.setGeometry(QtCore.QRect(260, 30, 50, 25))

        self.l_r2.setFont(font2)
        self.l_r2.setObjectName("l_r2")
        self.i_r2 = QtWidgets.QLineEdit(self.Clump)
        self.i_r2.setGeometry(QtCore.QRect(310, 30, 131, 25))

        self.i_r2.setFont(font1)
        self.i_r2.setObjectName("i_r2")
        self.l_kb = QtWidgets.QLabel(self.Clump)
        self.l_kb.setGeometry(QtCore.QRect(260, 60, 50, 25))

        self.l_kb.setFont(font2)
        self.l_kb.setObjectName("l_kb")
        self.i_kb = QtWidgets.QLineEdit(self.Clump)
        self.i_kb.setGeometry(QtCore.QRect(310, 60, 131, 25))

        self.i_kb.setFont(font1)
        self.i_kb.setObjectName("i_kb")
        self.Clump.raise_()
        self.PRS.raise_()
        self.b_default_prs.raise_()
        self.Function.addTab(self.QRS, "")
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
        font_t = font2
        font_t.setKerning(True)
        font_t.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.c_vep.setFont(font_t)
        self.c_vep.setObjectName("c_vep")
        self.PCA = QtWidgets.QGroupBox(self.c_sample_qc)
        self.PCA.setEnabled(True)
        self.PCA.setGeometry(QtCore.QRect(10, 130, 501, 141))
        self.PCA.setFont(font)
        self.PCA.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.PCA.setObjectName("PCA")
        self.l_pca_ref = QtWidgets.QLabel(self.PCA)
        self.l_pca_ref.setGeometry(QtCore.QRect(30, 25, 201, 25))

        self.l_pca_ref.setFont(font2)
        self.l_pca_ref.setObjectName("l_pca_ref")
        self.i_pca_ref = QtWidgets.QLineEdit(self.PCA)
        self.i_pca_ref.setGeometry(QtCore.QRect(164, 26, 251, 25))

        self.i_pca_ref.setFont(font1)
        self.i_pca_ref.setObjectName("i_pca_ref")
        self.s_pca_ref = QtWidgets.QToolButton(self.PCA)
        self.s_pca_ref.setGeometry(QtCore.QRect(410, 25, 41, 27))

        self.s_pca_ref.setFont(font1)
        self.s_pca_ref.setObjectName("s_pca_ref")
        self.l_pca_pop = QtWidgets.QLabel(self.PCA)
        self.l_pca_pop.setGeometry(QtCore.QRect(30, 65, 181, 25))

        self.l_pca_pop.setFont(font2)
        self.l_pca_pop.setObjectName("l_pca_pop")
        self.i_pop = QtWidgets.QLineEdit(self.PCA)
        self.i_pop.setGeometry(QtCore.QRect(370, 105, 81, 25))

        self.i_pop.setFont(font1)
        self.i_pop.setObjectName("i_pop")
        self.l_pop = QtWidgets.QLabel(self.PCA)
        self.l_pop.setGeometry(QtCore.QRect(269, 105, 101, 25))

        self.l_pop.setFont(font2)
        self.l_pop.setObjectName("l_pop")
        self.l_pop_sep = QtWidgets.QLabel(self.PCA)
        self.l_pop_sep.setGeometry(QtCore.QRect(30, 105, 61, 25))

        self.l_pop_sep.setFont(font2)
        self.l_pop_sep.setObjectName("l_pop_sep")
        self.s_pca_pop = QtWidgets.QToolButton(self.PCA)
        self.s_pca_pop.setGeometry(QtCore.QRect(410, 65, 41, 27))

        self.s_pca_pop.setFont(font1)
        self.s_pca_pop.setObjectName("s_pca_pop")
        self.i_pca_pop = QtWidgets.QLineEdit(self.PCA)
        self.i_pca_pop.setGeometry(QtCore.QRect(204, 66, 211, 25))

        self.i_pca_pop.setFont(font1)
        self.i_pca_pop.setObjectName("i_pca_pop")
        self.i_pop_id = QtWidgets.QLineEdit(self.PCA)
        self.i_pop_id.setGeometry(QtCore.QRect(210, 105, 41, 25))

        self.i_pop_id.setFont(font1)
        self.i_pop_id.setObjectName("i_pop_id")
        self.l_pop_id = QtWidgets.QLabel(self.PCA)
        self.l_pop_id.setGeometry(QtCore.QRect(160, 105, 41, 25))

        self.l_pop_id.setFont(font2)
        self.l_pop_id.setObjectName("l_pop_id")
        self.s_pop_sep = QtWidgets.QComboBox(self.PCA)
        self.s_pop_sep.setGeometry(QtCore.QRect(90, 105, 51, 25))

        self.s_pop_sep.setFont(font2)
        self.s_pop_sep.setObjectName("s_pop_sep")
        self.s_pop_sep.addItem("")
        self.s_pop_sep.addItem("")
        self.s_pop_sep.addItem("")
        self.l_pca_ref.raise_()
        self.i_pca_ref.raise_()
        self.s_pca_ref.raise_()
        self.l_pca_pop.raise_()
        self.i_pop.raise_()
        self.l_pop.raise_()
        self.l_pop_sep.raise_()
        self.i_pca_pop.raise_()
        self.i_pop_id.raise_()
        self.l_pop_id.raise_()
        self.s_pop_sep.raise_()
        self.s_pca_pop.raise_()
        self.MAF = QtWidgets.QGroupBox(self.c_sample_qc)
        self.MAF.setEnabled(True)
        self.MAF.setGeometry(QtCore.QRect(10, 60, 501, 61))
        self.MAF.setFont(font)
        self.MAF.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.MAF.setObjectName("MAF")
        self.l_maf_sample = QtWidgets.QLabel(self.MAF)
        self.l_maf_sample.setGeometry(QtCore.QRect(30, 25, 201, 25))

        self.l_maf_sample.setFont(font2)
        self.l_maf_sample.setObjectName("l_maf_sample")
        self.i_maf_sample = QtWidgets.QLineEdit(self.MAF)
        self.i_maf_sample.setGeometry(QtCore.QRect(164, 26, 251, 25))

        self.i_maf_sample.setFont(font1)
        self.i_maf_sample.setObjectName("i_maf_sample")
        self.s_maf_sample = QtWidgets.QToolButton(self.MAF)
        self.s_maf_sample.setGeometry(QtCore.QRect(410, 25, 41, 27))

        self.s_maf_sample.setFont(font1)
        self.s_maf_sample.setObjectName("s_maf_sample")
        self.Concordance = QtWidgets.QGroupBox(self.c_sample_qc)
        self.Concordance.setEnabled(True)
        self.Concordance.setGeometry(QtCore.QRect(10, 280, 501, 61))
        self.Concordance.setFont(font)
        self.Concordance.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.Concordance.setObjectName("Concordance")
        self.l_concord_ref = QtWidgets.QLabel(self.Concordance)
        self.l_concord_ref.setGeometry(QtCore.QRect(30, 25, 201, 25))

        self.l_concord_ref.setFont(font2)
        self.l_concord_ref.setObjectName("l_concord_ref")
        self.i_concord_ref = QtWidgets.QLineEdit(self.Concordance)
        self.i_concord_ref.setGeometry(QtCore.QRect(184, 26, 231, 25))

        self.i_concord_ref.setFont(font1)
        self.i_concord_ref.setObjectName("i_concord_ref")
        self.s_concord_ref = QtWidgets.QToolButton(self.Concordance)
        self.s_concord_ref.setGeometry(QtCore.QRect(410, 25, 41, 27))

        self.s_concord_ref.setFont(font1)
        self.s_concord_ref.setObjectName("s_concord_ref")
        self.b_default_sample_qc = QtWidgets.QPushButton(self.QC)
        self.b_default_sample_qc.setGeometry(QtCore.QRect(230, 380, 91, 31))
        self.b_default_sample_qc.setFont(font)
        self.b_default_sample_qc.setObjectName("b_default_sample_qc")
        self.Function.addTab(self.QC, "")
        self.Query_database = QtWidgets.QWidget()
        self.Query_database.setFont(font)
        self.Query_database.setObjectName("Query_database")
        self.c_query_database = QtWidgets.QGroupBox(self.Query_database)
        self.c_query_database.setGeometry(QtCore.QRect(30, 40, 501, 181))

        self.c_query_database.setFont(font)
        self.c_query_database.setCheckable(True)
        self.c_query_database.setObjectName("c_query_database")
        self.c_clinvar = QtWidgets.QCheckBox(self.c_query_database)
        self.c_clinvar.setGeometry(QtCore.QRect(40, 60, 451, 19))
        font_t = font
        font_t.setKerning(True)
        font_t.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.c_clinvar.setFont(font_t)
        self.c_clinvar.setChecked(True)
        self.c_clinvar.setObjectName("c_clinvar")
        self.c_pharmgkb = QtWidgets.QCheckBox(self.c_query_database)
        self.c_pharmgkb.setGeometry(QtCore.QRect(40, 110, 451, 19))

        self.c_pharmgkb.setFont(font_t)
        self.c_pharmgkb.setChecked(True)
        self.c_pharmgkb.setObjectName("c_pharmgkb")
        self.Function.addTab(self.Query_database, "")
        self.QR_code = QtWidgets.QWidget()

        self.QR_code.setFont(font)
        self.QR_code.setObjectName("QR_code")
        self.c_qr_code = QtWidgets.QGroupBox(self.QR_code)
        self.c_qr_code.setGeometry(QtCore.QRect(30, 40, 501, 341))

        self.c_qr_code.setFont(font)
        self.c_qr_code.setCheckable(True)
        self.c_qr_code.setObjectName("c_qr_code")
        self.l_qr_snps = QtWidgets.QLabel(self.c_qr_code)
        self.l_qr_snps.setGeometry(QtCore.QRect(30, 220, 281, 25))

        self.l_qr_snps.setFont(font)
        self.l_qr_snps.setObjectName("l_qr_snps")
        self.i_qr_snps = QtWidgets.QLineEdit(self.c_qr_code)
        self.i_qr_snps.setGeometry(QtCore.QRect(50, 261, 365, 25))

        self.i_qr_snps.setFont(font1)
        self.i_qr_snps.setObjectName("i_qr_snps")
        self.s_qr_snps = QtWidgets.QToolButton(self.c_qr_code)
        self.s_qr_snps.setGeometry(QtCore.QRect(410, 260, 41, 27))

        self.s_qr_snps.setFont(font1)
        self.s_qr_snps.setObjectName("s_qr_snps")
        self.l_or = QtWidgets.QLabel(self.c_qr_code)
        self.l_or.setGeometry(QtCore.QRect(230, 190, 31, 25))

        self.l_or.setFont(font)
        self.l_or.setObjectName("l_or")
        self.l_snps_list = QtWidgets.QLabel(self.c_qr_code)
        self.l_snps_list.setGeometry(QtCore.QRect(30, 30, 281, 25))

        self.l_snps_list.setFont(font)
        self.l_snps_list.setObjectName("l_snps_list")
        self.i_snps_list = QtWidgets.QPlainTextEdit(self.c_qr_code)
        self.i_snps_list.setGeometry(QtCore.QRect(30, 60, 421, 121))

        self.i_snps_list.setFont(font2)
        self.i_snps_list.setObjectName("i_snps_list")
        self.Function.addTab(self.QR_code, "")
        self.Advanced = QtWidgets.QWidget()
        self.Advanced.setObjectName("Advanced")
        self.c_ref_qc = QtWidgets.QGroupBox(self.Advanced)
        self.c_ref_qc.setGeometry(QtCore.QRect(20, 30, 521, 371))

        self.c_ref_qc.setFont(font)
        self.c_ref_qc.setCheckable(True)
        self.c_ref_qc.setObjectName("c_ref_qc")
        self.c_use_qc_ref = QtWidgets.QCheckBox(self.c_ref_qc)
        self.c_use_qc_ref.setGeometry(QtCore.QRect(30, 40, 451, 19))
        font_t = font2
        font_t.setKerning(True)
        font_t.setStyleStrategy(QtGui.QFont.PreferDefault)
        self.c_use_qc_ref.setFont(font_t)
        self.c_use_qc_ref.setChecked(False)
        self.c_use_qc_ref.setObjectName("c_use_qc_ref")
        self.Ref_qc_thresholds = QtWidgets.QGroupBox(self.c_ref_qc)
        self.Ref_qc_thresholds.setGeometry(QtCore.QRect(10, 80, 501, 271))
        font_t = QtGui.QFont()
        font_t.setFamily("Calibri")
        font_t.setBold(True)
        font_t.setWeight(75)
        self.Ref_qc_thresholds.setFont(font_t)
        self.Ref_qc_thresholds.setObjectName("Ref_qc_thresholds")
        self.l_qc_female = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_qc_female.setGeometry(QtCore.QRect(350, 110, 81, 25))

        self.l_qc_female.setFont(font2)
        self.l_qc_female.setToolTip("")
        self.l_qc_female.setObjectName("l_qc_female")
        self.l_vmiss = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_vmiss.setGeometry(QtCore.QRect(350, 70, 81, 25))

        self.l_vmiss.setFont(font2)
        self.l_vmiss.setToolTip("")
        self.l_vmiss.setObjectName("l_vmiss")
        self.l_het_mean = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_het_mean.setGeometry(QtCore.QRect(190, 190, 61, 25))

        self.l_het_mean.setFont(font2)
        self.l_het_mean.setObjectName("l_het_mean")
        self.l_hwe = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_hwe.setGeometry(QtCore.QRect(10, 150, 131, 25))

        self.l_hwe.setFont(font2)
        self.l_hwe.setToolTip("")
        self.l_hwe.setObjectName("l_hwe")
        self.i_hwe_p = QtWidgets.QLineEdit(self.Ref_qc_thresholds)
        self.i_hwe_p.setGeometry(QtCore.QRect(240, 150, 51, 25))

        self.i_hwe_p.setFont(font1)
        self.i_hwe_p.setObjectName("i_hwe_p")
        self.l_pihat = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_pihat.setGeometry(QtCore.QRect(10, 230, 161, 25))

        self.l_pihat.setFont(font2)
        self.l_pihat.setToolTip("")
        self.l_pihat.setObjectName("l_pihat")
        self.l_het = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_het.setGeometry(QtCore.QRect(10, 190, 151, 25))

        self.l_het.setFont(font2)
        self.l_het.setToolTip("")
        self.l_het.setObjectName("l_het")
        self.i_het = QtWidgets.QLineEdit(self.Ref_qc_thresholds)
        self.i_het.setGeometry(QtCore.QRect(260, 190, 31, 25))

        self.i_het.setFont(font1)
        self.i_het.setObjectName("i_het")
        self.i_f_female = QtWidgets.QLineEdit(self.Ref_qc_thresholds)
        self.i_f_female.setGeometry(QtCore.QRect(440, 110, 51, 25))

        self.i_f_female.setFont(font1)
        self.i_f_female.setObjectName("i_f_female")
        self.l_maf = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_maf.setGeometry(QtCore.QRect(10, 30, 171, 25))

        self.l_maf.setFont(font2)
        self.l_maf.setToolTip("")
        self.l_maf.setObjectName("l_maf")
        self.i_pihat = QtWidgets.QLineEdit(self.Ref_qc_thresholds)
        self.i_pihat.setGeometry(QtCore.QRect(190, 230, 141, 25))

        self.i_pihat.setFont(font1)
        self.i_pihat.setObjectName("i_pihat")
        self.i_smiss = QtWidgets.QLineEdit(self.Ref_qc_thresholds)
        self.i_smiss.setGeometry(QtCore.QRect(280, 70, 51, 25))

        self.i_smiss.setFont(font1)
        self.i_smiss.setObjectName("i_smiss")
        self.l_qc_male = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_qc_male.setGeometry(QtCore.QRect(190, 110, 81, 25))

        self.l_qc_male.setFont(font2)
        self.l_qc_male.setToolTip("")
        self.l_qc_male.setObjectName("l_qc_male")
        self.l_smiss = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_smiss.setGeometry(QtCore.QRect(190, 70, 81, 25))

        self.l_smiss.setFont(font2)
        self.l_smiss.setToolTip("")
        self.l_smiss.setObjectName("l_smiss")
        self.b_default_ref_qc = QtWidgets.QPushButton(self.Ref_qc_thresholds)
        self.b_default_ref_qc.setGeometry(QtCore.QRect(380, 200, 91, 31))

        self.b_default_ref_qc.setFont(font)
        self.b_default_ref_qc.setObjectName("b_default_ref_qc")
        self.i_f_male = QtWidgets.QLineEdit(self.Ref_qc_thresholds)
        self.i_f_male.setGeometry(QtCore.QRect(280, 110, 51, 25))

        self.i_f_male.setFont(font1)
        self.i_f_male.setObjectName("i_f_male")
        self.l_hwe_2 = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_hwe_2.setGeometry(QtCore.QRect(190, 150, 61, 25))

        self.l_hwe_2.setFont(font2)
        self.l_hwe_2.setObjectName("l_hwe_2")
        self.l_het_sd = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_het_sd.setGeometry(QtCore.QRect(300, 190, 41, 25))

        self.l_het_sd.setFont(font2)
        self.l_het_sd.setObjectName("l_het_sd")
        self.l_miss = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_miss.setGeometry(QtCore.QRect(10, 70, 181, 25))

        self.l_miss.setFont(font2)
        self.l_miss.setToolTip("")
        self.l_miss.setObjectName("l_miss")
        self.i_maf = QtWidgets.QLineEdit(self.Ref_qc_thresholds)
        self.i_maf.setGeometry(QtCore.QRect(190, 30, 141, 25))

        self.i_maf.setFont(font1)
        self.i_maf.setAccessibleDescription("")
        self.i_maf.setObjectName("i_maf")
        self.i_vmiss = QtWidgets.QLineEdit(self.Ref_qc_thresholds)
        self.i_vmiss.setGeometry(QtCore.QRect(440, 70, 51, 25))

        self.i_vmiss.setFont(font1)
        self.i_vmiss.setObjectName("i_vmiss")
        self.l_sex = QtWidgets.QLabel(self.Ref_qc_thresholds)
        self.l_sex.setGeometry(QtCore.QRect(10, 110, 181, 25))

        self.l_sex.setFont(font2)
        self.l_sex.setToolTip("")
        self.l_sex.setObjectName("l_sex")
        self.Ref_qc_thresholds.raise_()
        self.c_use_qc_ref.raise_()
        self.s_pca_pop.raise_()
        self.i_pop.raise_()
        self.Function.addTab(self.Advanced, "")
        PAGEANT.setCentralWidget(self.MainWindow)

        self.retranslateUi(PAGEANT)
        self.Function.setCurrentIndex(0)
        self.s_pop_sep.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(PAGEANT)

    def retranslateUi(self, PAGEANT):
        _translate = QtCore.QCoreApplication.translate
        PAGEANT.setWindowTitle(_translate("PAGEANT", "PAGEANT (2021-04-18)"))
        self.i_vcf.setToolTip(_translate("PAGEANT", "Select the genotype file that you want to analyze"))
        self.i_vcf.setText(_translate("PAGEANT", "./personal_genome/sample.vcf.gz"))
        self.b_analyze.setToolTip(_translate("PAGEANT", "Start analyze"))
        self.b_analyze.setText(_translate("PAGEANT", "Analyze"))
        self.s_vcf.setText(_translate("PAGEANT", "..."))
        self.l_vcf.setToolTip(_translate("PAGEANT", "Input the genotype file that you want to analyze"))
        self.l_vcf.setText(_translate("PAGEANT", "Genotype file"))
        self.l_output.setText(_translate("PAGEANT", "Ouput Directory"))
        self.i_output.setToolTip(_translate("PAGEANT", "Specify the output directory"))
        self.i_output.setText(_translate("PAGEANT", os.path.join(os.getcwd(), 'output')))
        self.s_output.setText(_translate("PAGEANT", "..."))
        self.l_name.setText(_translate("PAGEANT", "Name"))
        self.i_name.setToolTip(_translate("PAGEANT", "Input your name here"))
        self.i_name.setText(_translate("PAGEANT", "test"))
        self.Function.setTabText(self.Function.indexOf(self.Basic), _translate("PAGEANT", "IO"))
        self.s_ref.setText(_translate("PAGEANT", "..."))
        self.i_ref.setToolTip(_translate("PAGEANT", "Specify direcrtory of reference population data"))
        self.i_ref.setText(_translate("PAGEANT", "./population_genome"))
        self.l_database.setText(_translate("PAGEANT", "Database directory"))
        self.i_database.setToolTip(_translate("PAGEANT", "Specify directory of database"))
        self.i_database.setText(_translate("PAGEANT", "./algorithm_database"))
        self.s_database.setText(_translate("PAGEANT", "..."))
        self.l_ref.setText(_translate("PAGEANT", "Reference data directory"))
        self.b_default_prs.setToolTip(_translate("PAGEANT", "Set default values"))
        self.b_default_prs.setText(_translate("PAGEANT", "Default"))
        self.PRS.setTitle(_translate("PAGEANT", "PRS"))
        self.l_p_threshold.setToolTip(_translate("PAGEANT", "P threshold"))
        self.l_p_threshold.setText(_translate("PAGEANT", "P threshold"))
        self.l_ref_2.setText(_translate("PAGEANT", "Reference structure data"))
        self.l_EA.setToolTip(_translate("PAGEANT", "EA column name"))
        self.l_EA.setText(_translate("PAGEANT", "EA"))
        self.i_ld_ref.setToolTip(_translate("PAGEANT", "Specify directory of reference structure data"))
        self.i_ld_ref.setText(_translate("PAGEANT", os.path.abspath("./personal_genome/hapmap3.vcf.gz")))
        self.l_BETA.setToolTip(_translate("PAGEANT", "BETA column name"))
        self.l_BETA.setText(_translate("PAGEANT", "BETA"))
        self.i_p_threshold.setToolTip(_translate("PAGEANT", "P threshold for PRS"))
        self.i_p_threshold.setText(_translate("PAGEANT", "1e-5"))
        self.i_p_threshold.setPlaceholderText(_translate("PAGEANT", "P threshold for PRS"))
        self.s_ld_ref_2.setText(_translate("PAGEANT", "..."))
        self.i_EA.setToolTip(_translate("PAGEANT", "EA column name"))
        self.i_EA.setText(_translate("PAGEANT", "EA"))
        self.i_EA.setPlaceholderText(_translate("PAGEANT", "EA column name"))
        self.i_SNP.setToolTip(_translate("PAGEANT", "SNP column name"))
        self.i_SNP.setText(_translate("PAGEANT", "SNP"))
        self.i_SNP.setPlaceholderText(_translate("PAGEANT", "SNP column name"))
        self.l_P.setToolTip(_translate("PAGEANT", "P column name"))
        self.l_P.setText(_translate("PAGEANT", "P"))
        self.i_BETA.setToolTip(_translate("PAGEANT", "BETA column name"))
        self.i_BETA.setText(_translate("PAGEANT", "BETA"))
        self.i_BETA.setPlaceholderText(_translate("PAGEANT", "BETA column name"))
        self.i_P.setToolTip(_translate("PAGEANT", "P column name"))
        self.i_P.setText(_translate("PAGEANT", "P"))
        self.i_P.setPlaceholderText(_translate("PAGEANT", "P column name"))
        self.l_SNP.setToolTip(_translate("PAGEANT", "SNP column name"))
        self.l_SNP.setText(_translate("PAGEANT", "SNP"))
        self.Clump.setTitle(_translate("PAGEANT", "Clump"))
        self.l_p1.setToolTip(_translate("PAGEANT", "Index variant p-value threshold"))
        self.l_p1.setText(_translate("PAGEANT", "p1"))
        self.l_p2.setToolTip(_translate("PAGEANT", "Sites within clump p-value threshold"))
        self.l_p2.setText(_translate("PAGEANT", "p2"))
        self.i_p1.setToolTip(_translate("PAGEANT", "Index variant p-value threshold"))
        self.i_p1.setText(_translate("PAGEANT", "1"))
        self.i_p2.setToolTip(_translate("PAGEANT", "Sites within clump p-value threshold"))
        self.i_p2.setText(_translate("PAGEANT", "1"))
        self.l_r2.setToolTip(_translate("PAGEANT", "r^2 threshold"))
        self.l_r2.setText(_translate("PAGEANT", "r2"))
        self.i_r2.setToolTip(_translate("PAGEANT",
                                        "<html><head/><body><p>r<span style=\" vertical-align:super;\">2</span> threshold</p></body></html>"))
        self.i_r2.setText(_translate("PAGEANT", "0.1"))
        self.l_kb.setToolTip(_translate("PAGEANT", "Clump kb radius"))
        self.l_kb.setText(_translate("PAGEANT", "kb"))
        self.i_kb.setToolTip(_translate("PAGEANT", "Clump kb radius"))
        self.i_kb.setText(_translate("PAGEANT", "250"))
        self.i_kb.setPlaceholderText(_translate("PAGEANT", "Clump kb radius"))
        self.Function.setTabText(self.Function.indexOf(self.QRS), _translate("PAGEANT", "Quan traits"))
        # self.c_sample_qc.setToolTip(_translate("PAGEANT", "Activate function of sample QC"))
        self.c_sample_qc.setTitle(_translate("PAGEANT", "Sample QC"))
        self.c_vep.setToolTip(_translate("PAGEANT", "Using VEP to inspect the variants\' type in sample data"))
        self.c_vep.setText(_translate("PAGEANT", "Use VEP (Slow with too many variants or bad network!)"))
        self.PCA.setTitle(_translate("PAGEANT", "PCA plot"))
        self.l_pca_ref.setText(_translate("PAGEANT", "PCA reference data"))
        self.i_pca_ref.setToolTip(_translate("PAGEANT", "Specify the reference data for PCA"))
        self.i_pca_ref.setText(_translate("PAGEANT", "./personal_genome/hapmap3.vcf.gz"))
        self.i_pca_ref.setPlaceholderText(_translate("PAGEANT", "Select principal component analysis reference data"))
        self.s_pca_ref.setText(_translate("PAGEANT", "..."))
        self.l_pca_pop.setText(_translate("PAGEANT", "Reference population data"))
        self.i_pop.setToolTip(_translate("PAGEANT", "Separator in GWAS data"))
        self.i_pop.setText(_translate("PAGEANT", "population"))
        self.i_pop.setPlaceholderText(_translate("PAGEANT", "Input the column\' name of population"))
        self.l_pop.setText(_translate("PAGEANT", "Population Col."))
        self.l_pop_sep.setText(_translate("PAGEANT", "Text Sep."))
        self.s_pca_pop.setText(_translate("PAGEANT", "..."))
        self.i_pca_pop.setToolTip(_translate("PAGEANT", "Specify the reference population data for PCA"))
        self.i_pca_pop.setText(_translate("PAGEANT", "./personal_genome/hapmap3_samples.txt"))
        self.i_pca_pop.setPlaceholderText(_translate("PAGEANT", "Select principal component analysis population data"))
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
        self.Function.setTabText(self.Function.indexOf(self.QC), _translate("PAGEANT", "QC"))
        # self.c_query_database.setToolTip(_translate("PAGEANT", "Activate function of querying database"))
        self.c_query_database.setTitle(_translate("PAGEANT", "Query database"))
        self.c_clinvar.setToolTip(_translate("PAGEANT", "Query Clinvar Database"))
        self.c_clinvar.setText(_translate("PAGEANT", "Query Clinvar Database"))
        self.c_pharmgkb.setToolTip(_translate("PAGEANT", "Query PharmGKB Database"))
        self.c_pharmgkb.setText(_translate("PAGEANT", "Query PharmGKB Database"))
        self.Function.setTabText(self.Function.indexOf(self.Query_database), _translate("PAGEANT", "Query DB"))
        # self.c_qr_code.setToolTip(_translate("PAGEANT", "Activate function of generating QR code"))
        self.c_qr_code.setTitle(_translate("PAGEANT", "Generate QR code"))
        self.l_qr_snps.setText(_translate("PAGEANT", "Need SNPs\' list text"))
        self.i_qr_snps.setToolTip(_translate("PAGEANT", "Specify the text which include the needed SNPs list"))
        self.i_qr_snps.setText(_translate("PAGEANT", "./personal_genome/fingerprint_snps.txt"))
        self.s_qr_snps.setText(_translate("PAGEANT", "..."))
        self.l_or.setText(_translate("PAGEANT", "OR"))
        self.l_snps_list.setText(_translate("PAGEANT", "Need SNPs\' list"))
        self.i_snps_list.setPlaceholderText(
            _translate("PAGEANT", "Enter the list of SNPs, separated by space or new line "))
        self.Function.setTabText(self.Function.indexOf(self.QR_code), _translate("PAGEANT", "QR code"))
        # self.c_ref_qc.setToolTip(_translate("PAGEANT", "Activate function of reference QC"))
        self.c_ref_qc.setTitle(_translate("PAGEANT", "Population QC (Only support single reference file now)"))
        self.c_use_qc_ref.setToolTip(_translate("PAGEANT", "Use reference population data after QC in further analyze"))
        self.c_use_qc_ref.setText(_translate("PAGEANT", "Use population data after QC "))
        self.Ref_qc_thresholds.setTitle(_translate("PAGEANT", "Threshold for reference QC"))
        self.l_qc_female.setText(_translate("PAGEANT", "Female"))
        self.l_vmiss.setText(_translate("PAGEANT", "Per variant"))
        self.l_het_mean.setToolTip(_translate("PAGEANT", "r^2 threshold"))
        self.l_het_mean.setText(_translate("PAGEANT", "Mean  ±"))
        self.l_hwe.setText(_translate("PAGEANT", "P value of HWE"))
        self.i_hwe_p.setToolTip(_translate("PAGEANT", "The power of HWE\'s p-value maximum"))
        self.i_hwe_p.setText(_translate("PAGEANT", "50"))
        self.l_pihat.setText(_translate("PAGEANT", "Pihat of relatedness"))
        self.l_het.setText(_translate("PAGEANT", "Heterozygosity rate"))
        self.i_het.setToolTip(_translate("PAGEANT", "Heterozygosity rate outlier "))
        self.i_het.setText(_translate("PAGEANT", "3"))
        self.i_f_female.setToolTip(_translate("PAGEANT", "F of female minimum"))
        self.i_f_female.setText(_translate("PAGEANT", "0.6"))
        self.l_maf.setText(_translate("PAGEANT", "Minor allele frequency"))
        self.i_pihat.setToolTip(_translate("PAGEANT", "Pihat of relatedness maximum"))
        self.i_pihat.setText(_translate("PAGEANT", "0.2"))
        self.i_pihat.setPlaceholderText(_translate("PAGEANT", "Clump kb radius"))
        self.i_smiss.setToolTip(_translate("PAGEANT", "Missing genotype rate per sample maximum"))
        self.i_smiss.setText(_translate("PAGEANT", "0.02"))
        self.l_qc_male.setText(_translate("PAGEANT", "Male"))
        self.l_smiss.setText(_translate("PAGEANT", "Per sample"))
        self.b_default_ref_qc.setToolTip(_translate("PAGEANT", "Set default values"))
        self.b_default_ref_qc.setText(_translate("PAGEANT", "Default"))
        self.i_f_male.setToolTip(_translate("PAGEANT", "F of male maximum"))
        self.i_f_male.setText(_translate("PAGEANT", "0.4"))
        self.l_hwe_2.setToolTip(_translate("PAGEANT", "r^2 threshold"))
        self.l_hwe_2.setText(_translate("PAGEANT", "10 ^  -"))
        self.l_het_sd.setToolTip(_translate("PAGEANT", "r^2 threshold"))
        self.l_het_sd.setText(_translate("PAGEANT", "sd"))
        self.l_miss.setText(_translate("PAGEANT", "Missing genotype rate:"))
        self.i_maf.setToolTip(_translate("PAGEANT", "Minor allele frequency minimum"))
        self.i_maf.setText(_translate("PAGEANT", "0.01"))
        self.i_vmiss.setToolTip(_translate("PAGEANT", "Missing genotype rate per variant maximum"))
        self.i_vmiss.setText(_translate("PAGEANT", "0.02"))
        self.l_sex.setText(_translate("PAGEANT", "F of sex discrepancy:"))
        self.Function.setTabText(self.Function.indexOf(self.Advanced), _translate("PAGEANT", "Advanced"))


class MyMainForm(QMainWindow, Ui_PAGEANT):
    def __init__(self, parent=None):
        super(MyMainForm, self).__init__(parent)
        self.setupUi(self)

        self.s_vcf.clicked.connect(self.o_vcf)
        self.s_output.clicked.connect(self.o_output)
        self.s_database.clicked.connect(self.o_database)
        self.s_ref.clicked.connect(self.o_ref)
        self.s_ld_ref_2.clicked.connect(self.o_ld_ref)
        self.s_qr_snps.clicked.connect(self.o_qr_list)
        self.s_maf_sample.clicked.connect(self.o_maf)
        self.s_pca_ref.clicked.connect(self.o_pca)
        self.s_pca_pop.clicked.connect(self.o_text)
        self.s_concord_ref.clicked.connect(self.o_concord)
        self.b_analyze.clicked.connect(self.run)
        self.b_default_prs.clicked.connect(self.set_default_prs)
        self.b_default_sample_qc.clicked.connect(self.set_default_sample_qc)
        self.b_default_ref_qc.clicked.connect(self.set_default_ref_qc)

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

    def o_output(self):
        get_dir_path = QFileDialog.getExistingDirectory(self, "Select the output directory", os.getcwd())
        self.i_output.setText(get_dir_path)
        self.parameters['output'] = get_dir_path

    def o_database(self):
        get_dir_path = QFileDialog.getExistingDirectory(self, "Select the database directory", os.getcwd())
        self.i_database.setText(get_dir_path)
        self.parameters['database'] = get_dir_path

    def o_ref(self):
        get_dir_path = QFileDialog.getExistingDirectory(self, "Select the reference population directory", os.getcwd())
        self.i_ref.setText(get_dir_path)
        self.parameters['ref'] = get_dir_path

    def o_ld_ref(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the reference structure file", os.getcwd(),
                                                         filter='vcf (*.vcf *.vcf.gz)')
        self.i_ld_ref.setText(get_directory_path[0])
        self.parameters['ref_structure'] = get_directory_path[0]

    def o_maf(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the MAF reference file", os.getcwd(),
                                                         filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; '
                                                                'PLINK 2 binary (*.pgen)')
        self.i_maf_sample.setText(get_directory_path[0])
        self.parameters['maf_ref'] = get_directory_path[0]

    def o_pca(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the PCA reference file", os.getcwd(),
                                                         filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; '
                                                                'PLINK 2 binary (*.pgen)')
        self.i_pca_ref.setText(get_directory_path[0])
        self.parameters['pca_ref'] = get_directory_path[0]

    def o_concord(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the concordance analysis reference file", os.getcwd(),
                                                         filter='vcf (*.vcf *.vcf.gz);; PLINK 1 binary (*.bed);; '
                                                                'PLINK 2 binary (*.pgen)')
        self.i_concord_ref.setText(get_directory_path[0])
        self.parameters['concord_ref'] = get_directory_path[0]

    def o_qr_list(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the SNPs list file", os.getcwd())
        self.i_qr_snps.setText(get_directory_path[0])
        self.parameters['qr_snps_txt'] = get_directory_path[0]

    def o_text(self):
        get_directory_path = QFileDialog.getOpenFileName(self, "Select the PCA population file", os.getcwd())
        self.i_pca_pop.setText(get_directory_path[0])
        self.parameters['population_file'] = get_directory_path[0]

    def set_default_prs(self):
        self.i_p1.setText(str(default_config['clump-p1']))
        self.i_p2.setText(str(default_config['clump-p2']))
        self.i_r2.setText(str(default_config['clump-r2']))
        self.i_kb.setText(str(default_config['clump-kb']))
        self.i_P.setText(default_config['P'])
        self.i_BETA.setText(default_config['BETA'])
        self.i_SNP.setText(default_config['SNP'])
        self.i_EA.setText(default_config['EA'])
        self.i_p_threshold.setText(str(default_config['p_threshold']))
        self.i_ld_ref.setText(default_config['ref_structure'])

    def set_default_ref_qc(self):
        self.i_maf.setText(str(default_config['qc_maf']))
        self.i_smiss.setText(str(default_config['qc_smiss']))
        self.i_vmiss.setText(str(default_config['qc_vmiss']))
        self.i_hwe_p.setText(str(default_config['qc_hardy']))
        self.i_pihat.setText(str(default_config['qc_pihat']))
        self.i_f_male.setText(str(default_config['qc_male_F']))
        self.i_f_female.setText(str(default_config['qc_female_F']))
        self.i_het.setText(str(default_config['qc_het']))

    def set_default_sample_qc(self):
        self.i_maf_sample.setText(default_config['maf_ref'])
        self.i_pca_ref.setText(str(default_config['pca_ref']))
        self.i_pca_pop.setText(default_config['population_file'])
        self.i_concord_ref.setText(str(default_config['concord_ref']))
        self.i_pop.setText(default_config['population_col'])
        self.i_pop_id.setText(str(default_config['population_id']))
        self.s_pop_sep.setCurrentIndex(0)

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
        self.parameters.update(dict(zip(['name', 'data_dir', 'input_file', 'ref', 'ref_structure', 'output',
                                         'SNP', 'EA', 'P', 'BETA', 'p_threshold', 'clump-p1', 'clump-r2',
                                         'clump-kb', 'clump-p2', 'qr_snps_txt', 'qc_maf', 'qc_vmiss', 'qc_smiss',
                                         'qc_hardy', 'qc_het', 'qc_male_F', 'qc_female_F', 'qc_pihat', 'use_qc_ref',
                                         'sample_qc', 'ref_qc', 'vep', 'query_database', 'pharmgkb',
                                         'clinvar', 'qr_code', 'maf_ref', 'pca_ref', 'concord_ref',
                                         'population_sep', 'population_col', 'population_file', 'population_id'
                                         ],
                                        [self.i_name.text(), self.i_database.text(), self.i_vcf.text(),
                                         self.i_ref.text(), self.i_ld_ref.text(), self.i_output.text(),
                                         self.i_SNP.text(), self.i_EA.text(), self.i_P.text(), self.i_BETA.text(),
                                         float(self.i_p_threshold.text()), float(self.i_p1.text()),
                                         float(self.i_r2.text()), int(self.i_kb.text()), float(self.i_p2.text()),
                                         self.i_qr_snps.text(), float(self.i_maf.text()), float(self.i_vmiss.text()),
                                         float(self.i_smiss.text()), int(self.i_hwe_p.text()), int(self.i_het.text()),
                                         float(self.i_f_male.text()), float(self.i_f_female.text()),
                                         float(self.i_pihat.text()), self.c_use_qc_ref.isChecked(),
                                         self.c_sample_qc.isChecked(), self.c_ref_qc.isChecked(),
                                         self.c_vep.isChecked(), self.c_query_database.isChecked(),
                                         self.c_pharmgkb.isChecked(), self.c_clinvar.isChecked(),
                                         self.c_qr_code.isChecked(),self.i_maf_sample.text(), self.i_pca_ref.text(),
                                         self.i_concord_ref.text(), index_to_sep(self.s_pop_sep.currentIndex()),
                                         self.i_pop.text(), self.i_pca_pop.text(), self.i_pop_id.text()
                                         ])))
        if self.i_snps_list.toPlainText():
            self.parameters.update({'qr_snps_txt': self.i_snps_list.toPlainText()})
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
