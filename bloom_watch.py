from qgis.PyQt.QtWidgets import QAction
from qgis.PyQt.QtGui import QIcon
from .ui_handler import UIHandler
import os


class BloomWatchPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.ui = None
        self.action = None

    def initGui(self):
        icon_path = os.path.join(os.path.dirname(__file__), "icons", "plugin_icon.png")
        self.action = QAction(
            QIcon(icon_path), "NASA Data Loader", self.iface.mainWindow()
        )
        self.action.triggered.connect(self.run)
        self.iface.addPluginToMenu("NASA", self.action)
        self.iface.addToolBarIcon(self.action)

    def unload(self):
        self.iface.removePluginMenu("NASA", self.action)
        self.iface.removeToolBarIcon(self.action)

    def run(self):
        self.ui = UIHandler(self.iface)
        self.ui.show()
        self.ui.activateWindow()
