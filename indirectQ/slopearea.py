# USGS methods for indirect discharge calculations

import pandas as pd
import geopandas as gpd


class Sections:

    def __init__(self, siteID, agency):

        self.siteID = siteID
        self.agency = agency
        # FUTURE DEVELOPMENT
        #   Define methods for looking up location information using input site ID and operating agency
        self.siteName = None
        self.siteLat = None
        self.siteLon = None
        # geometry variables
        self.XSections = None
        self.Centerline = None
        self.HWMs = None

    def load_Xsections(self, in_paths):


    def load_ManningsRaster(self, man_path):


    def load_topography(self, topo_path):

    def extract_Xsections(self):

    def _

    def calculateQ(self):