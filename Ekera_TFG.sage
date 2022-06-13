#! /usr/bin/sage -python
# coding: utf-8
# Usar factorizar_completamente(N, c = 1, k = None).

# Clase auxiliar para guardar los factores no triviales de N en forma reducida.
class ColeccionFactores:
    def __init__(self, N):
        # El entero N a factorizar
        self.N = N

        # El conjunto de factores encontrado hasta ahora, reducido para que todos los factores del conjunto
        # sean coprimos por pares entre sí. Esta propiedad se aplica mediante add().
        self.factores_encontrados = set()

        # El conjunto de factores primos encontrados hasta ahora; un subconjunto de factores_encontrados
        self.primos_encontrados = set()

        # El residuo; el producto de los factores coprimos por pares compuestos en el
        # colección, o uno si no hay factores compuestos en la colección.
        self.residuo = 1

        # Añade N como factor.
        self.add(N)

    # Comprueba si todos los factores primos se han encontrado.
    def esta_completa(self):
        return self.residuo == 1

    # Añade un factor a la colección.
    def add(self, d):
        # Verificamos que el factor no sea trivial y que aún no se haya encontrado.
        if (d == 1) or (d in self.factores_encontrados):
            return

        # Comprobamos si d comparte un factor con cualquiera de los factores encontrados.
        D = 1

        for f in self.factores_encontrados:
            D = gcd(f, d)

            if D != 1:
                break

        if D != 1:
            # Si es así, eliminamos f, dividimos f y d, y añadimos los factores resultantes.
            self.factores_encontrados.remove(f)
            if f not in self.primos_encontrados:
                # También eliminamos f del residuo cuando eliminamos f de la colección.
                self.residuo /= f

            f /= D
            d /= D

            self.add(D);

            if f != 1:
                self.add(f)

            if d != 1:
                self.add(d);
        else:
            # Comprobamos si d es una potencia perfecta y, de ser así, reducimos d.
            (d, _) = ZZ(d).perfect_power();

            # Añadimos d a los factores encontrados
            self.factores_encontrados.add(d);

            # Comprobamos si d es primo, y si es así lo añadimos
            resultado = d.is_prime(proof = False);

            if resultado:
                self.primos_encontrados.add(d);
            else:
                # Si d no es primo, multiplicamos d por el residuo.
                self.residuo *= d;

    # Imprimimos la información de los conjuntos
    def imprimir_estado(self):
        print("Factores encontrados:", len(self.factores_encontrados));
        print("Pirmos encontrados:", len(self.primos_encontrados));

        factores_encontrados = list(self.factores_encontrados);
        factores_encontrados.sort();

        for i in range(len(factores_encontrados)):
            print(" Factor " + str(i) + ":", factores_encontrados[i]);
        print("");


# ------------------------------------------------------------------------------
# Función que factoriza completamente N (nos da todos sus factores primos sin multiplicidad)
#
# El parámetro c esta definido en el trabajo, es una constante mas grande o igual que 1. El parámetro k
# no necesita especificarse explícitamente: por defecto, tantas iteraciones k como sean necesarias para
# factorizar completamente el entero N.
#
# Si se quiere, se puede especificar k. Si el número de
# iteraciones realizadas excede k, entonces se detiene el algoritmo
#
# Esta función devuelve el conjunto de todos los factores primos distintos que dividen a N.

def factorizar_completamente(N, c = 1, k = None):

    # Comprobaciones de los parametros
    if (N < 2) or (c < 1):
        raise Exception("Error: Parameteros incorrectos.")

    # Función de soporte para construir el producto de q^e, para q primos <= B y
    # e el exponente mas grande tal que q^e <= B para una cota B.
    def productorio_potencias_primas(B):
        factor = 1
        for q in prime_range(B + 1):
            e = 1
            while q^(e + 1) <= B:
                e += 1
            factor *= q^e

        return factor

    # Función de soporte para calcular t tal que x = 2^t * o para o impar.
    def exponente_t(x):
        if x == 0:
            return 0

        t = 0
        while (x % 2) == 0:
            t += 1
            x /= 2

        return t

    #Función auxiliar que dado N, calcula un g aleatoriamente de Z_N^* y nos da su orden
    def orden(N):
        while True:
            g = IntegerModRing(N).random_element() #Seleccionamos g aleatoriamente de Z_N^*
            if (g == 1):
                continue
            if gcd(g.lift(), N) == 1:
                break
        return g.multiplicative_order() #Devolvemos el orden de g (Esto es lo que nos daria Shor)

    # Paso 1: Seleccionamos g aleatoriamente de Z_N^* y calculamos su orden
    r = orden(N)

    # Paso 2: Construimos el producto de factores primos q^e < cm y lo multiplicamos por r.
    # Obtenemos r' que denominamos rp

    m = N.nbits() # m es el número de bits de N
    rp = productorio_potencias_primas(c * m) * r #Calculamos r'

    # Paso 3: Sea rp = 2^t o para o impar.
    t = exponente_t(rp) #calculamos t
    o = rp / 2^t #obtenemos o

    # Definimos el conjunto de coprimos a pares y añadimos N.
    F = ColeccionFactores(N);

    # Paso 4: Para j = 1, 2, ... hasta k donde k puede estar o no acotada.
    j = 0;
    while True:
        # Imprimimos el estado actual para cada iteración
        print("Iteracion:", j);
        F.imprimir_estado();

        # Comprobamos si ya hemos acabado
        if F.esta_completa():
            break;

        # Incrementamos j para la siguiente iteración.
        j += 1;

        # Comprobamos que j > k, si k se ha especificado, y en ese caso
        # devolvemos una excepción.
        if (k != None) and (j > k):
            raise Exception("Error: Se ha superado el límite de iteraciones.");

        # Paso 4.1: Seleccionamos x uniformemente al azar de Z_N^*.
        x=0
        while x == 0:
            x = IntegerModRing(N).random_element();
        x = IntegerModRing(N)(x) #Para aplicar exponenciacion/aritmetica modular

        #4.2 Para cada i = 0, ..., t hacemos:
        for i in range(0, t + 1):
            xj = x^((2^i) * o);
            # Paso 4.2.1 para i = 0, 1, ..., t:
            d = gcd((xj - 1).lift(), N);
            if 1 < d < N:
                F.add(d);

    return F.primos_encontrados;
