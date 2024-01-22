
import numpy as np
from PyQt5 import QtWidgets
import sys
import scipy
from scipy import special
import math
import numpy as np
import matplotlib.pyplot as plt

class Widget():

    def __init__(self):
        #Блок с основными параметрами задающими приложение
        app = QtWidgets.QApplication(sys.argv)
        self.window = QtWidgets.QWidget()
        self.window1 = QtWidgets.QWidget()
        self.window2 = QtWidgets.QWidget()
        self.window.resize(600,400)
        self.window1.resize(400,400)
        self.window2.resize(400, 400)
        #=======================================#
        #Настройка первого окна
        #=======================================#
        # Настройка закладок
        self.tabb = QtWidgets.QTabWidget()
        self.tabb.addTab(self.window1, "Программа один ")
        self.tabb.addTab(self.window2, "Программа два")
        # =======================================#
        #Блок с массивами параметров графиков для распределения Парето
        self.paramsPareto = []
        self.paramsParetoNear1 = []
        self.paramsParetoNear0 = []
        self.legendPareto = []
        self.legendParetoNear1 = []
        self.legendParetoNear0 = []
        # =======================================#
        #Блок с массивами параметров графиков для распределения Гнеденко-Вейбула
        self.paramsGnedenkoVeybula = []
        self.paramsGnedenkoVeybulaNear1 = []
        self.paramsGnedenkoVeybulaNear0 = []
        self.legendGnedenkoVeybula = []
        self.legendGnedenkoVeybulaNear1 = []
        self.legendGnedenkoVeybulaNear0 = []
        # =======================================#
        # Блок с массивами параметров графиков для Логнормального распределения
        self.paramsLogNorm = []
        self.paramsLogNormNear1 = []
        self.paramsLogNormNear0 = []
        self.legendLogNorm = []
        self.legendLogNormNear1 = []
        self.legendLogNormNear0 = []
        #=======================================#
        #Конект кнопок с их слотами (которое потом меняется в зависимости от выбранного распределения)
        self.btnDraw1 = QtWidgets.QPushButton("Draw Pareto")
        self.btnDraw2 = QtWidgets.QPushButton("Draw Pareto near 1")
        self.btnDraw3 = QtWidgets.QPushButton("Draw Pareto near 0")
        self.btnClearALl = QtWidgets.QPushButton("Clear All")
        self.btnDraw1.clicked.connect(self.slotBtnPareto)
        self.btnDraw2.clicked.connect(self.slotBtnParetoNear1)
        self.btnDraw3.clicked.connect(self.slotBtnParetoNear0)
        self.btnClearALl.clicked.connect(self.clearAll)
        # =======================================#
        # Настройка ComboBox, добавление туда элементов
        self.combobox = QtWidgets.QComboBox()
        self.indexComboBox = 0
        self.combobox.currentIndexChanged.connect(self.slotComboBox)
        self.combobox.insertItem(0,"Распределение Парето")
        self.combobox.insertItem(1,"Распределение Гнеденко-Вейбула")
        self.combobox.insertItem(2, "Логнормальное распределение")

        #=======================================#
        #Настройка внешних параметров и окна ввода
        self.lineV = QtWidgets.QLineEdit()
        labelV = QtWidgets.QLabel("Коэффициент вариации v = ")

        hboxV = QtWidgets.QHBoxLayout()

        hboxV.addWidget(labelV)
        hboxV.addWidget(self.lineV)

        # =======================================#
        #Добавление всего в layout и настройка приложения под layout
        vbox = QtWidgets.QVBoxLayout()

        vbox.addLayout(hboxV)
        vbox.addWidget(self.combobox)
        vbox.addWidget(self.btnDraw1)
        vbox.addWidget(self.btnDraw2)
        vbox.addWidget(self.btnDraw3)
        vbox.addWidget(self.btnClearALl)

        self.window1.setLayout(vbox)
        #=======================================#
        #Настройка окна для второй программы

        # =======================================#
        # Настройка ComboBox, добавление туда элементов Window 2

        self.w2_comboboxA = QtWidgets.QComboBox()
        self.w2_indexComboBoxA = 0
        self.w2_comboboxA.currentIndexChanged.connect(self.w2_slotComboBoxA)
        self.w2_comboboxA.insertItem(0, "Распределение Гнеденко-Вейбулла")
        self.w2_comboboxA.insertItem(1, "Логнормальное распределение")
        self.w2_comboboxA.insertItem(2, "Гамма распределение")

        self.w2_comboboxB = QtWidgets.QComboBox()
        self.w2_indexComboBoxB = 0
        self.w2_comboboxB.currentIndexChanged.connect(self.w2_slotComboBoxB)
        self.w2_comboboxB.insertItem(0, "Распределение Гнеденко-Вейбулла")
        self.w2_comboboxB.insertItem(1, "Логнормальное распределение")
        self.w2_comboboxB.insertItem(2, "Гамма распределение")


        # =======================================#
        # Настройка вывода и ввода Window 2

        self.JIJA = QtWidgets.QLineEdit()

        self.w2_labelPValue = QtWidgets.QLabel("----PValue----")

        self.w2_labelEpsValue = QtWidgets.QLabel("----Eps-Value----")
        self.w2_lineEps = QtWidgets.QLineEdit()
        self.w2_lineEps.setText("0.0001")

        self.w2_lineMA = QtWidgets.QLineEdit()
        self.w2_lineMA.setText("1")
        w2_labelMA = QtWidgets.QLabel("Математического ожидание первого распределения mu1 = ")
        self.w2_lineMB = QtWidgets.QLineEdit()
        self.w2_lineMB.setText("0.01")
        w2_labelMB = QtWidgets.QLabel("Математическое ожидание второго распределения mu2 = ")
        self.w2_lineVA = QtWidgets.QLineEdit()
        self.w2_lineVA.setText("0.5")
        w2_labelVA = QtWidgets.QLabel("Коэффициент вариации первого распределения v1 = ")
        self.w2_lineVB = QtWidgets.QLineEdit()
        self.w2_lineVB.setText("0.5")
        w2_labelVB = QtWidgets.QLabel("Коэффициент вариации второго распределения v2 = ")

        w2_layout = QtWidgets.QVBoxLayout()

        w2_hboxEps = QtWidgets.QHBoxLayout()

        w2_btnLayout = QtWidgets.QHBoxLayout()
        w2_hboxAM = QtWidgets.QHBoxLayout()
        w2_hboxAV = QtWidgets.QHBoxLayout()
        w2_vboxA = QtWidgets.QVBoxLayout()

        w2_hboxBM = QtWidgets.QHBoxLayout()
        w2_hboxBV = QtWidgets.QHBoxLayout()
        w2_vboxB = QtWidgets.QVBoxLayout()
        # =======================================#
        # Конект кнопок с их слотами Window 2

        self.w2_btnGetPValue = QtWidgets.QPushButton("Расчитать вероятность")
        self.w2_btnGetPValue.clicked.connect(self.w2_slotBtn)
        # =======================================#
        # Layout для распределения А
        w2_hboxAM.addWidget(w2_labelMA)
        w2_hboxAM.addWidget(self.w2_lineMA)
        w2_hboxAV.addWidget(w2_labelVA)
        w2_hboxAV.addWidget(self.w2_lineVA)
        w2_vboxA.addLayout(w2_hboxAM)
        w2_vboxA.addLayout(w2_hboxAV)
        w2_vboxA.addWidget(self.w2_comboboxA)
        # =======================================#
        # Layout для распределения В
        w2_hboxBM.addWidget(w2_labelMB)
        w2_hboxBM.addWidget(self.w2_lineMB)
        w2_hboxBV.addWidget(w2_labelVB)
        w2_hboxBV.addWidget(self.w2_lineVB)
        w2_vboxB.addLayout(w2_hboxBM)
        w2_vboxB.addLayout(w2_hboxBV)
        w2_vboxB.addWidget(self.w2_comboboxB)
        # =======================================#
        # Основной Layout
        w2_layout.addLayout(w2_vboxA)
        w2_layout.addLayout(w2_vboxB)
        # =======================================#
        # Layout для кнопки и точности
        w2_btnLayout.addWidget(self.w2_btnGetPValue)
        w2_btnLayout.addWidget(self.w2_labelPValue)
        w2_hboxEps.addWidget(self.w2_labelEpsValue)
        w2_hboxEps.addWidget(self.w2_lineEps)
        # =======================================#
        #Основной Layout
        w2_layout.addLayout(w2_hboxEps)
        w2_layout.addLayout(w2_btnLayout)


        self.window2.setLayout(w2_layout)

        #=======================================#
        # Настройка закрытия приложения
        self.btnQuit = QtWidgets.QPushButton("Exit")
        self.btnQuit.clicked.connect(app.quit)
        # =======================================#
        # Добавление закладок в основной layout Window 2
        main_vbox = QtWidgets.QVBoxLayout()
        main_vbox.addWidget(self.tabb)
        main_vbox.addWidget(self.btnQuit)
        # =======================================#
        # Настройка основного окна
        self.window.setLayout(main_vbox)
        #Открыте/Закрытие окна
        self.window.show()
        sys.exit(app.exec_())


        # =======================================#
        # =======================================#


    # =======================================#
    # Блок отвечающий за ComboBox и смену работы кнопок Window 1#
    # =======================================#
    def slotComboBox(self, index):
        self.clearAll()
        self.btnDraw1.clicked.disconnect()
        self.btnDraw2.clicked.disconnect()
        self.btnDraw3.clicked.disconnect()
        if (index == 0):
            self.btnDraw1.setText("Draw Pareto")
            self.btnDraw2.setText("Draw Pareto near 1")
            self.btnDraw3.setText("Draw Pareto near 0")
            self.btnDraw1.clicked.connect(self.slotBtnPareto)
            self.btnDraw2.clicked.connect(self.slotBtnParetoNear1)
            self.btnDraw3.clicked.connect(self.slotBtnParetoNear0)
        elif(index == 1):
            self.btnDraw1.setText("Draw Gnedenko-Veybula")
            self.btnDraw2.setText("Draw Gnedenko-Veybula near 1")
            self.btnDraw3.setText("Draw Gnedenko-Veybula near 0")
            self.btnDraw1.clicked.connect(self.slotBtnGnedenkoVeybula)
            self.btnDraw2.clicked.connect(self.slotBtnGnedenkoVeybulaNear1)
            self.btnDraw3.clicked.connect(self.slotBtnGnedenkoVeybulaNear0)
        elif(index == 2):
            self.btnDraw1.setText("Draw LogNorm")
            self.btnDraw2.setText("Draw LogNorm near 1")
            self.btnDraw3.setText("Draw LogNorm near 0")
            self.btnDraw1.clicked.connect(self.slotBtnLogNorm)
            self.btnDraw2.clicked.connect(self.slotBtnLogNormNear1)
            self.btnDraw3.clicked.connect(self.slotBtnLogNormNear0)

    # =======================================#
    # Конец блока отвечающего за ComboBox и смену работы кнопок#
    # =======================================#



    # =======================================#
    # Блок с распределением Парето#
    # =======================================#

    # =======================================#
    # Функция отклонения безотказной работы R(tv) для распределения Парето зависящая от коэффициента вариации#
    def fPareto(self,t,v):
        bf = math.sqrt(1 + 1/v**2)
        return (bf/(t*(1+bf)))**(1+bf)

    # =======================================#
    # Так как распределение Парето имеет отступ, этот отступ надо найти , тогда можно получить аргумент и саму фунцию#
    def getDrawParamsPareto(self,v):
        up = math.sqrt(1.0 + 1 / v ** 2)
        t_0 = up / (1.0 + up) # Начальное время распределения Парето
        t = np.linspace(t_0, 10, 1000)
        graph = [self.fPareto(i, v) for i in t]
        return [t,graph]

    # =======================================#
    # Можно удалить#
    def drawPareto(self):
        v = float(self.lineV.text())
        print(self.getDrawParamsPareto(v))

    # =======================================#
    # Слот отвечающий за отрисовку графиков распределения Парето с введенным коэффициентом вариации#
    def slotBtnPareto(self):
        v = float(self.lineV.text())#
        [t,Ftv] = self.getDrawParamsPareto(v)
        labelGraph = "c коэфф. вариации v = " + str(v)

        location = ['upper right', 'upper left', 'lower left',#локации на графике ,где могут находится легенды
                    'lower right', 'right', 'center left',
                    'center right', 'lower center', 'upper center', 'center']
        self.legendPareto.append(labelGraph)#Добавление подписи для легенд
        self.paramsPareto.append([t, Ftv])#и параметров для графика
        plt.ioff() ## Позволяет убрать предыдущий график и отчертить все заново
        plt.close()## на новом графике (добавив новые функции)
        fig, ax = plt.subplots()
        for i in range(len(self.paramsPareto)):
            tt = self.paramsPareto[i][0]#Отрисовка всех графиков ранее добавленных в наши массивы
            Ftt = self.paramsPareto[i][1]
            ax.plot(tt,Ftt, label = self.legendPareto[i])
            ax.legend(loc = location[0])
            print(self.legendPareto[i])
        # plt.plot(t,Ftv,'-g')
        plt.title("Функция больших уклонений для распределения Парето")
        plt.show()
        self.combobox.showPopup()

    def slotBtnParetoNear1(self):
        v = float(self.lineV.text())
        [t,Ftv] = self.getDrawParamsPareto(v)
        nt = []
        nFtv = []
        for i in range(20):
            nt.append(t[i])
            nFtv.append(Ftv[i])
        t = nt
        Ftv = nFtv
        labelGraph = "c коэфф. вариации v = " + str(v)
        location = ['upper right', 'upper left', 'lower left',
                    'lower right', 'right', 'center left',
                    'center right', 'lower center', 'upper center', 'center']
        self.legendParetoNear1.append(labelGraph)
        self.paramsParetoNear1.append([t, Ftv])
        plt.ioff() ## Позволяет убрать предыдущий график и отчертить все заново
        plt.close()## на новом графике (добавив новые функции)
        fig, ax = plt.subplots()
        for i in range(len(self.paramsParetoNear1)):
            tt = self.paramsParetoNear1[i][0]
            Ftt = self.paramsParetoNear1[i][1]
            ax.plot(tt,Ftt, label = self.legendParetoNear1[i])
            ax.legend(loc = location[2])
        # plt.plot(t,Ftv,'-g')
        plt.title("Функция больших уклонений для распределeния Парето")
        plt.show()

    def slotBtnParetoNear0(self):
        v = float(self.lineV.text())
        [t,Ftv] = self.getDrawParamsPareto(v)
        nt = []
        nFtv = []
        for i in range(100,200):
            nt.append(t[i])
            nFtv.append(Ftv[i])
        t = nt
        Ftv = nFtv
        labelGraph = "c коэфф. вариации v = " + str(v)
        location = ['upper right', 'upper left', 'lower left',
                    'lower right', 'right', 'center left',
                    'center right', 'lower center', 'upper center', 'center']
        self.legendParetoNear0.append(labelGraph)
        self.paramsParetoNear0.append([t, Ftv])
        plt.ioff() ## Позволяет убрать предыдущий график и отчертить все заново
        plt.close()## на новом графике (добавив новые функции)
        fig, ax = plt.subplots()
        for i in range(len(self.paramsParetoNear0)):
            tt = self.paramsParetoNear0[i][0]
            Ftt = self.paramsParetoNear0[i][1]
            ax.plot(tt,Ftt, label = self.legendParetoNear0[i])
            ax.legend(loc = location[0])
        # plt.plot(t,Ftv,'-g')
        plt.title("Функция больших уклонений для распределeния Парето")
        plt.show()


    # =======================================#
    # Конец блока с распределением Парето#
    # =======================================#

    # =======================================#
    # Блок с распределением Гнеденко-Вейбула#
    # =======================================#
    def fGnedenkoVeybula(self,t,a):
        return (math.exp(-(t)**a*(scipy.special.gamma(1+1/a))**a))

    def getDrawGnedenkoVeybula(self,v):
        a = self.solveGnedenkoVeybula(v)
        t = np.linspace(0, 10, 1000)
        graph = [self.fGnedenkoVeybula(i, a) for i in t]
        return [t,graph]

    def slotBtnGnedenkoVeybula(self):
        v = float(self.lineV.text())
        [t , Ftv] = self.getDrawGnedenkoVeybula(v)
        labelGraph = "c коэфф. вариации v = " + str(v)
        location = ['upper right', 'upper left', 'lower left',
                    'lower right', 'right', 'center left',
                    'center right', 'lower center', 'upper center', 'center']
        self.legendGnedenkoVeybula.append(labelGraph)  # Добавление подписи для легенд
        self.paramsGnedenkoVeybula.append([t, Ftv])  # и параметров для графика
        plt.ioff()  ## Позволяет убрать предыдущий график и отчертить все заново
        plt.close()  ## на новом графике (добавив новые функции)
        fig, ax = plt.subplots()
        for i in range(len(self.paramsGnedenkoVeybula)):
            tt = self.paramsGnedenkoVeybula[i][0]  # Отрисовка всех графиков ранее добавленных в наши массивы
            Ftt = self.paramsGnedenkoVeybula[i][1]
            ax.plot(tt, Ftt, label=self.legendGnedenkoVeybula[i])
            ax.legend(loc=location[0])
        # plt.plot(t,Ftv,'-g')
        plt.title("Функция больших уклонений для распределения Гнеденко-Вейбула")
        plt.show()

    def slotBtnGnedenkoVeybulaNear1(self):
        v = float(self.lineV.text())
        [t , Ftv] = self.getDrawGnedenkoVeybula(v)
        nt = []
        nFtv = []
        for i in range(10):
            nt.append(t[i])
            nFtv.append(Ftv[i])
        t = nt
        Ftv = nFtv
        labelGraph = "c коэфф. вариации v = " + str(v)
        location = ['upper right', 'upper left', 'lower left',
                    'lower right', 'right', 'center left',
                    'center right', 'lower center', 'upper center', 'center']
        self.legendGnedenkoVeybulaNear1.append(labelGraph)  # Добавление подписи для легенд
        self.paramsGnedenkoVeybulaNear1.append([t, Ftv])  # и параметров для графика
        plt.ioff()  ## Позволяет убрать предыдущий график и отчертить все заново
        plt.close()  ## на новом графике (добавив новые функции)
        fig, ax = plt.subplots()
        for i in range(len(self.paramsGnedenkoVeybulaNear1)):
            tt = self.paramsGnedenkoVeybulaNear1[i][0]  # Отрисовка всех графиков ранее добавленных в наши массивы
            Ftt = self.paramsGnedenkoVeybulaNear1[i][1]
            ax.plot(tt, Ftt, label=self.legendGnedenkoVeybulaNear1[i])
            ax.legend(loc=location[0])
            print(self.legendGnedenkoVeybulaNear1[i])
        # plt.plot(t,Ftv,'-g')
        plt.title("Функция больших уклонений для распределения Гнеденко-Вейбула")
        plt.show()

    def slotBtnGnedenkoVeybulaNear0(self):
        v = float(self.lineV.text())
        [t , Ftv] = self.getDrawGnedenkoVeybula(v)
        nt = []
        nFtv = []
        for i in range(100,200):
            nt.append(t[i])
            nFtv.append(Ftv[i])
        t = nt
        Ftv = nFtv
        labelGraph = "c коэфф. вариации v = " + str(v)
        location = ['upper right', 'upper left', 'lower left',
                    'lower right', 'right', 'center left',
                    'center right', 'lower center', 'upper center', 'center']
        self.legendGnedenkoVeybulaNear0.append(labelGraph)  # Добавление подписи для легенд
        self.paramsGnedenkoVeybulaNear0.append([t, Ftv])  # и параметров для графика
        plt.ioff()  ## Позволяет убрать предыдущий график и отчертить все заново
        plt.close()  ## на новом графике (добавив новые функции)
        fig, ax = plt.subplots()
        for i in range(len(self.paramsGnedenkoVeybulaNear0)):
            tt = self.paramsGnedenkoVeybulaNear0[i][0]  # Отрисовка всех графиков ранее добавленных в наши массивы
            Ftt = self.paramsGnedenkoVeybulaNear0[i][1]
            ax.plot(tt, Ftt, label=self.legendGnedenkoVeybulaNear0[i])
            ax.legend(loc=location[0])
            print(self.legendGnedenkoVeybulaNear0[i])
        # plt.plot(t,Ftv,'-g')
        plt.title("Функция больших уклонений для распределения Гнеденко-Вейбула")
        plt.show()

    #========Вспомогательные функции=========#
    def vGnedenkoVeybula(self,a):
        return (scipy.special.gamma(1 + 2 / a) / (scipy.special.gamma(1 + 1 / a)) ** 2 - 1) ** (1. / 2)

    def solveGnedenkoVeybula(self,v):
        a1 = 0.015
        a2 = 10.0
        while (math.fabs(self.vGnedenkoVeybula(a1) - v) > 0.0001):
            if (self.vGnedenkoVeybula((a1 + a2) / 2) - v > 0):
                a1 = (a1 + a2) / 2
            else:
                a2 = (a1 + a2) / 2
        print("Параметр В = " + str(a1))
        return(a1)

    # =======================================#
    # Конец блока с распределением Гнеденко-Вейбула#
    # =======================================#

    # =======================================#
    # Блок с Логнормальным рапределением#
    # =======================================#
    def fLogNorm(self, t, v):
        up = math.log(t*v*math.exp(1/2))
        down = math.sqrt(2*math.log(v**2 + 1))

        return 1/2 - (1/2)*special.erf(up/down)

    def getDrawLogNorm(self,v):
        t0 = 1/(v*math.sqrt(math.exp(1/2)))#Начальное время функции отклонения
        t = np.linspace(0.0001, 10, 1000)
        graph = [self.fLogNorm(i, v) for i in t]
        return [t,graph]

    def slotBtnLogNorm(self):
        v = float(self.lineV.text())
        [t , Ftv] = self.getDrawLogNorm(v)
        labelGraph = "c коэфф. вариации v = " + str(v)
        location = ['upper right', 'upper left', 'lower left',
                    'lower right', 'right', 'center left',
                    'center right', 'lower center', 'upper center', 'center']
        self.legendLogNorm.append(labelGraph)  # Добавление подписи для легенд
        self.paramsLogNorm.append([t, Ftv])  # и параметров для графика
        plt.ioff()  ## Позволяет убрать предыдущий график и отчертить все заново
        plt.close()  ## на новом графике (добавив новые функции)
        fig, ax = plt.subplots()
        for i in range(len(self.paramsLogNorm)):
            tt = self.paramsLogNorm[i][0]  # Отрисовка всех графиков ранее добавленных в наши массивы
            Ftt = self.paramsLogNorm[i][1]
            ax.plot(tt, Ftt, label=self.legendLogNorm[i])
            ax.legend(loc=location[0])
        # plt.plot(t,Ftv,'-g')
        plt.title("Функция больших уклонений для Логнормального распределения")
        plt.show()

    def slotBtnLogNormNear1(self):
        v = float(self.lineV.text())
        [t , Ftv] = self.getDrawLogNorm(v)
        labelGraph = "c коэфф. вариации v = " + str(v)
        location = ['upper right', 'upper left', 'lower left',
                    'lower right', 'right', 'center left',
                    'center right', 'lower center', 'upper center', 'center']
        nt = []
        nFtv = []
        for i in range(10):
            nt.append(t[i])
            nFtv.append(Ftv[i])
        t = nt
        Ftv = nFtv
        self.legendLogNormNear1.append(labelGraph)  # Добавление подписи для легенд
        self.paramsLogNormNear1.append([t, Ftv])  # и параметров для графика
        plt.ioff()  ## Позволяет убрать предыдущий график и отчертить все заново
        plt.close()  ## на новом графике (добавив новые функции)
        fig, ax = plt.subplots()
        for i in range(len(self.paramsLogNormNear1)):
            tt = self.paramsLogNormNear1[i][0]  # Отрисовка всех графиков ранее добавленных в наши массивы
            Ftt = self.paramsLogNormNear1[i][1]
            ax.plot(tt, Ftt, label=self.legendLogNormNear1[i])
            ax.legend(loc=location[0])
        # plt.plot(t,Ftv,'-g')
        plt.title("Функция больших уклонений для Логнормального распределения")
        plt.show()

    def slotBtnLogNormNear0(self):
        v = float(self.lineV.text())
        [t , Ftv] = self.getDrawLogNorm(v)
        labelGraph = "c коэфф. вариации v = " + str(v)
        location = ['upper right', 'upper left', 'lower left',
                    'lower right', 'right', 'center left',
                    'center right', 'lower center', 'upper center', 'center']
        nt = []
        nFtv = []
        for i in range(200,400):
            nt.append(t[i])
            nFtv.append(Ftv[i])
        t = nt
        Ftv = nFtv
        self.legendLogNormNear0.append(labelGraph)  # Добавление подписи для легенд
        self.paramsLogNormNear0.append([t, Ftv])  # и параметров для графика
        plt.ioff()  ## Позволяет убрать предыдущий график и отчертить все заново
        plt.close()  ## на новом графике (добавив новые функции)
        fig, ax = plt.subplots()
        for i in range(len(self.paramsLogNormNear0)):
            tt = self.paramsLogNormNear0[i][0]  # Отрисовка всех графиков ранее добавленных в наши массивы
            Ftt = self.paramsLogNormNear0[i][1]
            ax.plot(tt, Ftt, label=self.legendLogNormNear0[i])
            ax.legend(loc=location[0])
        # plt.plot(t,Ftv,'-g')
        plt.title("Функция больших уклонений для Логнормального распределения")
        plt.show()

    # =======================================#
    # Конец блока с Логнормальным рапределением Window 1#
    # =======================================#

    # =======================================#
    # Блок отвечающий за очистку массивов Window 1#
    # =======================================#
    def clearAll(self):
        self.paramsPareto.clear()
        self.paramsParetoNear1.clear()
        self.paramsParetoNear0.clear()
        self.legendPareto.clear()
        self.legendParetoNear1.clear()
        self.legendParetoNear0.clear()

        self.paramsGnedenkoVeybula.clear()
        self.paramsGnedenkoVeybulaNear1.clear()
        self.paramsGnedenkoVeybulaNear0.clear()
        self.legendGnedenkoVeybula.clear()
        self.legendGnedenkoVeybulaNear1.clear()
        self.legendGnedenkoVeybulaNear0.clear()

        self.paramsLogNorm.clear()
        self.paramsLogNormNear1.clear()
        self.paramsLogNormNear0.clear()
        self.legendLogNorm.clear()
        self.legendLogNormNear1.clear()
        self.legendLogNormNear0.clear()

    # =======================================#
    # Блок отвечающий за ComboBox и смену работы кнопок Window 2#
    # =======================================#
    def w2_slotComboBoxA(self, index):
        self.w2_indexComboBoxA = index
        print("Индекс распределения А - " + str(self.w2_indexComboBoxA))

    def w2_slotComboBoxB(self, index):
        self.w2_indexComboBoxB = index
        print("Индекс распределения В - " + str(self.w2_indexComboBoxB))

    def w2_slotBtn(self):
        MA = float(self.w2_lineMA.text())
        VA = float(self.w2_lineVA.text())
        IndA = self.w2_indexComboBoxA
        IndB = self.w2_indexComboBoxB
        MB = float(self.w2_lineMB.text())
        VB = float(self.w2_lineVB.text())
        Eps = float(self.w2_lineEps.text())
        if(MA*VA*MB*VB*Eps == 0):
            print("Не все параметры заполнены")
            #Здесь можно добавить всплывающее окно
            self.w2_labelPValue.setText("----PValue----")
            return 0
        P = self.getPvalue(MA,VA,MB,VB,IndA,IndB, Eps)
        labelPV =  "----" + str(P) + "----"
        self.w2_labelPValue.setText(labelPV)


    # =======================================#
    # Блок с поиском интеграла (вероятности быть в состоянии 2 отказа)#
    # =======================================#
    def getPvalue(self , A_mu , A_v , B_mu , B_v , A_ind , B_ind , eps):
        if (A_ind > 3 | B_ind > 3):
            print("Головой подумай")
        A_ind = A_ind + 1
        B_ind = B_ind + 1
        #Распределение для А
        x1 = 0.0
        Acdf = lambda x,A,B: 1 - math.exp(-A*x**B)
        if (A_ind == 1):
            #Распределение Гнеденко-Вейбулла
            A_par_2 = self.solveGnedenkoVeybula(A_v)
            A_par_1 = ((scipy.special.gamma(1 + 1/A_par_2))/A_mu)**A_par_2
        elif (A_ind == 2):
            # Распределение Логнормальное распределение
            A_par_2 = math.sqrt(A_v**2 + 1);
            A_par_1 = math.log(A_mu) - (A_par_2**2)/2;

            Acdf = lambda x, mu, sigma: 1/2 + 1/2*special.erf((math.log(x) - mu) / (math.sqrt(2)*sigma))
        elif (A_ind == 3):
            # Здесь будет Гамма распреледелние
            Acdf = lambda x, A, B: 1 - math.exp(-A * x ** B)
        #Плотность для В
        if (B_ind == 1):
            #Распределение Гнеденко-Вейбулла
            B_par_2 = self.solveGnedenkoVeybula(B_v)
            B_par_1 = ((scipy.special.gamma(1 + 1/B_par_2))/B_mu)**B_par_2
            print("B1 = " + str(B_par_1) + "| B2 = " + str(B_par_2))
            #Квантиль нужного уровня
            x1 = ((math.log(1/(1 - eps/2)))/B_par_1)**(1/B_par_2)
            x_buff = ((math.log(1/(0.95)))/B_par_1)**(1/B_par_2)
            x2 = ((math.log(1/(eps/2)))/B_par_1)**(1/B_par_2)
            print("eps = " + str(eps) +"| x1 = " + str(x1) + "| x2 = " + str(x2))
            Bpdf = lambda x,A,B: (A*B*x**(B-1))*math.exp(-A*x**B)

        elif (B_ind == 2):
            # Распределение Логнормальное распределение
            B_par_2 = math.sqrt(B_v**2 + 1);
            B_par_1 = math.log(B_mu) - (B_par_2**2)/2;
            print("B1 = " + str(B_par_1) + "| B2 = " + str(B_par_2))
            #Квантиль нужного уровня
            x1 = 0.000000001
            x1_b = math.exp(math.sqrt(2)*B_par_2*scipy.special.erfinv(eps/2) + B_par_1)
            x_buff = math.exp(math.sqrt(2)*B_par_2*scipy.special.erfinv(0.05) + B_par_1)

            if (x1_b < x1):
                print("Произошла замена")
                x1 = x1_b
            x2 = math.exp(math.sqrt(2)*B_par_2*scipy.special.erfinv(1-eps) + B_par_1)
            # print("eps = " + str(eps) +"| x1 = " + str(x1) + "| x2 = " + str(x2))
            Bpdf = lambda x, mu, sigma:  math.exp(-((math.log(x) - mu) / (math.sqrt(2)*sigma))**2)/(math.sqrt(2*math.pi)*x*sigma)
        elif (B_ind == 3):
            # Здесь будет Гамма распреледелние
            Bpdf = lambda x, A, B: 1 - math.exp(-A * x ** B)
        L = round(math.fabs(x2 - x1)) + 1
        # добавить разбиение сетки!!!
        # добавить разбиение сетки!!!
        # добавить разбиение сетки!!! возле нуля
        # добавить разбиение сетки!!!
        # добавить разбиение сетки!!!

        n1 = 100000
        # n1 = round((1/eps)) + 1
        n2 = L*100
        # x_buff = x1 + (x2 - x1)/100
        xx1 = np.linspace(x1 , x_buff ,n1)
        xx2 = np.linspace(x_buff,x2,n2)
        x = np.hstack((xx1,xx2))
        n = n1 + n2
        x = np.linspace(x1,x2,n)
        f = lambda x: Acdf(x,A_par_1,A_par_2)*Bpdf(x,B_par_1,B_par_2)
        g = lambda x: Bpdf(x,B_par_1,B_par_2)
        q = 0
        Intq = 0
        for i in range(n - 1):
            q = q + (x[i+1] - x[i])*(f(x[i]) + 3*f((2*x[i] + x[i+1])/3) + 3*f((x[i] + 2*x[i+1])/3) + f(x[i+1]))/8
            Intq = Intq + (x[i+1] - x[i])*(g(x[i]) + 3*g((2*x[i] + x[i+1])/3) + 3*g((x[i] + 2*x[i+1])/3) + g(x[i+1]))/8
            # if (x[i+1] - x[i] < 0):
                # print("x[i+1] | " + str(x[i+1]) + " || x[i] | " + str(x[i]))
        C_mu = B_mu*3
        print("____________________________________")
        print("q | " + str(q))
        print("Значение интеграла| " + str(Intq))
        pi_2 = C_mu*q/(A_mu + q*(A_mu+C_mu))
        print("pi_2 | " + str(pi_2))
        return pi_2
    # =======================================#
    # Конец блока с поиском интеграла (вероятности быть в состоянии 2 отказа)#
    # =======================================#



if __name__ == "__main__":
    widget = Widget()
    # print(math.log(math.exp(1)))
    # print(math.exp(1))
    # t = [1,1]
    # tt = [[1,1], [2,2]]
    # ttt = [tt , [[1,2],[1,2]],[[1,2],[1,2]]]
    # print(all([t in jija for jija in ttt]))