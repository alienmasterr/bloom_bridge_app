from PyQt5.QtWidgets import (
    QAbstractItemView, QPushButton, QMessageBox, QHBoxLayout, QDialog, QWidget
)
from PyQt5.QtCore import QTimer, QDateTime, Qt, QSettings
from qgis.core import QgsProject
from PyQt5.QtGui import QFont, QIcon
from PyQt5 import uic
import os


try:
    from alina_processor import *
except Exception as e:
    print(f"❌ Failed to import alina_processor: {e}")
    main = None

FORM_CLASS, _ = uic.loadUiType(os.path.join(os.path.dirname(__file__), "form.ui"))

class UIHandler(QDialog, FORM_CLASS):
    def __init__(self, iface):
        super().__init__()
        self.iface = iface
        self.setupUi(self)
        self.setFont(QFont("Segoe UI", 10))
        font = QFont("Segoe UI", 10)

        self.dataset_select.setFont(font)
        self.location_select.setFont(font)
        self.date_from.setFont(font)
        self.date_to.setFont(font)
        self.new_layer_name.setFont(font)
        self.existing_layer_select.setFont(font)
        self.progress_status.setFont(font)
        self.progress_bar.setFont(font)

        self.setWindowIcon(QIcon(os.path.join(os.path.dirname(__file__), "icons", "plugin_icon.png")))
        self.setStyleSheet("QDialog { background-color: #fff; }")
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setFormat("%p%")

        self.dataset_select.setSelectionMode(QAbstractItemView.MultiSelection)
        self.location_select.setSelectionMode(QAbstractItemView.MultiSelection)

        self.dataset_select.itemSelectionChanged.connect(self.on_dataset_selected)
        self.radio_new_layer.toggled.connect(self.update_layer_fields)
        self.radio_existing_layer.toggled.connect(self.update_layer_fields)
        self.load_button.clicked.connect(self.on_load_clicked)

        self.load_session()
        self.populate_datasets()
        self.populate_existing_layers()

    def populate_datasets(self):
        self.dataset_select.clear()
        for layer in QgsProject.instance().mapLayers().values():
            if layer.type() == layer.VectorLayer:
                self.dataset_select.addItem(layer.name())

    def populate_existing_layers(self):
        self.existing_layer_select.clear()
        for layer in QgsProject.instance().mapLayers().values():
            if layer.type() == layer.RasterLayer:
                self.existing_layer_select.addItem(layer.name())

    def update_layer_fields(self):
        self.new_layer_name.setEnabled(self.radio_new_layer.isChecked())
        self.existing_layer_select.setEnabled(self.radio_existing_layer.isChecked())

    def on_dataset_selected(self):
        selected_datasets = [item.text() for item in self.dataset_select.selectedItems()]
        if not selected_datasets:
            return

        self.date_from.setEnabled(True)
        self.date_to.setEnabled(True)
        self.location_select.setEnabled(True)
        self.load_target_radio.setEnabled(True)
        self.radio_new_layer.setEnabled(True)
        self.radio_existing_layer.setEnabled(True)

        dates, locations = self.extract_metadata(selected_datasets)
        if dates:
            self.date_from.setDateTime(min(dates))
            self.date_to.setDateTime(max(dates))

        self.location_select.clear()
        for loc in locations:
            self.location_select.addItem(loc)

        default_name = f"{selected_datasets[0]}_{self.date_from.dateTime().toString('dd.MM.yyyy')}-{self.date_to.dateTime().toString('dd.MM.yyyy')}"
        self.new_layer_name.setText(default_name)

    def on_load_clicked(self):
        locations = [item.text() for item in self.location_select.selectedItems()]
        if not locations:
            self.progress_status.setText("⚠️ No location selected")
            return

        layer_name = (
            self.new_layer_name.text()
            if self.radio_new_layer.isChecked()
            else self.existing_layer_select.currentText()
        )

        self.progress_status.setText("🚀 Processing started...")
        self.load_button.setEnabled(False)
        self.progress_bar.setMaximum(0)
        self.progress_bar.setValue(0)

        try:
            if main:
                main()  # ← Alina's processing script
                self.progress_status.setText("✅ Processing complete")
            else:
                raise ImportError("alina_processor module not imported")
        except Exception as e:
            self.progress_status.setText(f"❌ Processing error: {str(e)}")

        self.load_button.setEnabled(True)
        self.save_session()

    def extract_metadata(self, layer_names):
        dates = []
        locations = set()
        for layer in QgsProject.instance().mapLayers().values():
            if layer.name() in layer_names:
                for f in layer.getFeatures():
                    date_val = f.attribute("date_x")
                    lat = f.attribute("decimalLatitude")
                    lon = f.attribute("decimalLongitude")
                    if date_val:
                        if isinstance(date_val, QDateTime):
                            dates.append(date_val)
                        elif isinstance(date_val, str):
                            dt = QDateTime.fromString(date_val, Qt.ISODate)
                            if dt.isValid():
                                dates.append(dt)
                    if lat and lon:
                        locations.add(f"{lat}, {lon}")
        return dates, sorted(locations)

    def save_session(self):
        settings = QSettings()
        settings.setValue("nasa_loader/datasets", [item.text() for item in self.dataset_select.selectedItems()])
        settings.setValue("nasa_loader/date_from", self.date_from.dateTime().toString(Qt.ISODate))
        settings.setValue("nasa_loader/date_to", self.date_to.dateTime().toString(Qt.ISODate))
        settings.setValue("nasa_loader/locations", [item.text() for item in self.location_select.selectedItems()])

    def load_session(self):
        settings = QSettings()
        datasets = settings.value("nasa_loader/datasets", [])
        date_from = settings.value("nasa_loader/date_from")
        date_to = settings.value("nasa_loader/date_to")
        locations = settings.value("nasa_loader/locations", [])

        for i in range(self.dataset_select.count()):
            item = self.dataset_select.item(i)
            if item.text() in datasets:
                item.setSelected(True)

        if date_from:
            self.date_from.setDateTime(QDateTime.fromString(str(date_from), Qt.ISODate))
        if date_to:
            self.date_to.setDateTime(QDateTime.fromString(str(date_to), Qt.ISODate))

        for i in range(self.location_select.count()):
            item = self.location_select.item(i)
            if item.text() in locations:
                item.setSelected(True)