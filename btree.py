class BTreeNode:

    #key id of the node
    #parent is a BTreeNode
    #children is a list of BTreeNodes
    def __init__(self, key, parent, children):
        self.__key = key
        self.__parent = parent
        self.__children = []

    def addChild(self, child):
        self.__children.append(child)
    
    def isLeaf(self):
        return len(self.__children) == 0
            
    def isRoot(self):
        return self.__parent == None

class BTree:

    #rootNode is a BTreeNode
    def __init__(self, rootNode):
        self.__root = rootNode
        self.__numNode = 1
        
