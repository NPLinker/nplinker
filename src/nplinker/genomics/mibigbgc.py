from .bgc import BGC


class MiBIGBGC(BGC):

    def __init__(self, id, strain, name, product_prediction):
        super().__init__(id, strain, name, product_prediction)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id

