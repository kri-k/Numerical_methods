import sys
from PySide import QtGui
from libs.widget_control import get_widget_from_ui
import lab_5_1.gui.main_window
# from lab_5_2.gui import *
# from lab_5_3.gui import *
# from lab_5_4.gui import *

# from libs.arithmetic_parser import FunctionFromStr


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    main_menu = get_widget_from_ui('./gui/ui_forms/main_menu.ui')

    lab_5_1_widget = lab_5_1.gui.main_window.MainWindow(parent=main_menu)
    main_menu.pushButton_parabolic.clicked.connect(lab_5_1_widget.showMaximized)

    main_menu.pushButton_hyperbolic.setEnabled(False)
    main_menu.pushButton_elliptical.setEnabled(False)
    main_menu.pushButton_parabolic2D.setEnabled(False)

    main_menu.show()
    sys.exit(app.exec_())
