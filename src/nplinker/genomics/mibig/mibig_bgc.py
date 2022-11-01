from ..bgc import BGC


class MibigBGC(BGC):

    def __init__(self, id: int, strain: str, name: str, product_prediction: str):
        """The object reprenseting MiBIG BGC annotations/metadata

        Args:
            id(int): used to identify a BGC in a group of BGCs
            strain(str): strain name
            name(str): BGC name
            product_prediction(str): prediction of the product encoded by BGC
        """

        super().__init__(id, strain, name, product_prediction)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id
