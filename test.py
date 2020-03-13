import scipy as sp

class Foo:

	def __init__(self):
		self.list = []

	def getK(self, num):
		self.list = sp.zeros(num)

foo = Foo()


print(foo.list, type(foo.list))

foo.getK(5)

print(foo.list, type(foo.list))



		

