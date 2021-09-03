'''
Clase: Entidad que define una serie de elementos que determinan un estado (atributos) y un comportamiento (método). Plantillas de los objetos.
Objeto: Instancia de una clase
'''

class animal():
    ## Use constructor
    def __init__(self, nombre, edad, ruido):
        self.nombre = nombre
        self.edad = edad
        self.ruido = ruido
    def haz_ruido(self):
        print(self.ruido)

perro = animal('Rex', 10, 'XOXO')
perro.haz_ruido()
print(perro.__dict__)

'''
Herencia: Mecanismo por el cual una clase se deriva de otra de manera que extiende su funcionalidad. Al derivarse de una clase, hereda sus métodos y atributos.
Superclase y subclase: La clase de la que estamos heredando se suele denominar clase base, clase padre, superclase, etc. La clase que recibe la herencia es la subclase.
'''

class perro(animal):
  pass

class gato(animal):
  usa_arenero = True

Aeden = gato('Aeden', 9, 'Miau')
print(Aeden.__dict__)
Rex = perro('Rex', 10, 'Guau')
print(Rex.__dict__)

'''
Overriding: Implementar la información heredada de forma diferente.
'''

class perro_guau(perro):
  def cambiar_ruido(self):
    self.ruido = 'guau'

Eskel = perro_guau('Eskel', 13, 'J')
print(Eskel.__dict__)
Eskel.cambiar_ruido()
print(Eskel.__dict__)

'''
Polimorfismos: Al invocar un método con el mismo nombre en distintos objetos obtendremos una respuesta distinta.
'''

class gato(gato):
  def haz_ruido(self):
      self.ruido = 'Miau'
      print(self.ruido)

Ves = gato('Ves', 14, 'KlKl')
print(Ves.__dict__)
Ves.haz_ruido()
print(Ves.__dict__)

'''
BioProjects: bioJulia, Bioconductor, BioRuby, BioPerl, BioJava
Bibliotecas:
NumPy: Arrays, multi-D, vectorization, Broadcasting, numeros imaginarios.
Pandas (Panel Data): Dataframes, indexing integrado, estructura de datos.
Matplotlib: Emplea Pyplot, lo que da interfaz similar a MATLAB, conectado a NumPy y Pandas, exploración de analisis de datos y plots para publicaciones, problemas con datasets grandes, visualización interactiva para web y graphical user interfaces.
Seaborn: Librería para hacer gráficas estadísticas en Python, construida en matplotlib, integra estructuras de datos de Pandas y NumPy.
SciPy: Colección de algoritmos matemáticos y funciones creadas en extensión de NumPy.
Fuentes:
http://biopython.org/DIST/docs/tutorial/Tutorial.html
https://numpy.org/doc/stable/user/whatisnumpy.html
https://pandas.pydata.org/about/index.html
https://matplotlib.org/stable/users/index.html
https://seaborn.pydata.org/introduction.html
https://docs.scipy.org/doc/scipy/reference/tutorial/general.html
'''
