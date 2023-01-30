from extractors import CrisprCorporateExtractor


class AGenomeExtractorFactory:

    def create_extractor(self, genome):
        pass


class GenomeExtractorFactory(AGenomeExtractorFactory):

    def create_extractor(self, genome):
        extractor = None
        if genome == "human":
            pass
        else:
            extractor = CrisprCorporateExtractor()

        return extractor
