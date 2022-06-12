# Algoritmo-de-Shor-y-variantes
TFG Juan Cano - UB 2022
Contenidos:
1. Notebook de Jupyter con la explicación y realización del tutorial de Qiskit de IBM, con comentarios de cada linia y de los resultados y alternativas de implementación. Kernel Python 3.8.3
2. Notebook de Jupyter con la realización del algoritmo descrito en la memoria de la variante de Ekera del algoritmo de Shor basado en la obra de Miller. Realizado con SageMath 9.2
3. Este README con la explicación del algoritmo de Shor y un ejemplo para el caso $N=21$ y $a=2$
Notebook 

### De la factorización a la determinación del período

La teoría de números que subyace en el algoritmo de Shor se relaciona con secuencias de módulo-periódicas. Echemos un vistazo a un ejemplo de tal secuencia. Consideremos la secuencia de las potencias de dos:
$$1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, ...$$
Ahora veamos la misma secuencia 'módulo 15', es decir, el resto después de dividir por quince cada una de estas potencias de dos:
$$1, 2, 4, 8, 1, 2, 4, 8, 1, 2, 4, ...$$
Esta es una secuencia módulo 15 que se repite cada cuatro números, es decir, es una secuencia módulo-periódica con un período de cuatro.

La reducción de la factorización de $N$ al problema de encontrar el período de un entero $x$ menor que $N$ y mayor que $1$ depende del siguiente resultado de la teoría de números:

> La función $\mathcal{F}(a) = x^a \bmod N$ es una función periódica, donde $x$ es un coprimo entero de $N$ y $a \ge 0$.

Hay que tener en cuenta que dos números son coprimos, si el único entero positivo que los divide a ambos es 1. Esto equivale a que su máximo común divisor sea 1. Por ejemplo, 8 y 15 son coprimos, ya que no comparten ningún factor común (distinto de 1). Sin embargo, 9 y 15 no son coprimos, ya que ambos son divisibles por 3.

> Dado que $\mathcal{F}(a)$ es una función periódica, por lo tanto, tiene algún período $r$. Sabiendo que $x^0 \bmod N = 1$, esto significa que $x^r \bmod N = 1$ ya que la función es periódica y, por lo tanto, $r$ es la primera potencia distinta de cero donde $x^r = 1 (\bmod N)$.

Dada esta información y mediante la siguiente manipulación algebraica:
$$ x^r \equiv 1 \bmod N $$
$$ x^r = (x^{r/2})^2 \equiv 1 \bmod N $$
$$ (x^{r/2})^2 - 1 \equiv 0 \bmod N $$
y si $r$ es un número par:
$$ (x^{r/2} + 1)(x^{r/2} - 1) \equiv 0 \bmod N $$

A partir de esto, el producto $(x^{r/2} + 1)(x^{r/2} - 1)$ es un múltiplo entero de $N$, el número a factorizar. Por lo tanto, siempre que $(x^{r/2} + 1)$ o $(x^{r/2} - 1)$ no sea un múltiplo de $N$, entonces al menos uno de $(x^ {r/2} + 1)$ o $(x^{r/2} - 1)$ debe tener un factor no trivial en común con $N$.

Por lo que al calcular $\text{gcd}(x^{r/2} - 1, N)$ y $\text{gcd}(x^{r/2} + 1, N)$ se obtendrá un factor de $N$, donde $\text{mcd}$ es la función del máximo común denominador, que se puede calcular en tiempo polinomial.

#### Pasos clásicos del algoritmo de Shor

Para simplificar, supongamos que $N$ tiene solo dos factores primos distintos: $N = pq$.

<div class="alert alert-block alert-info"> <a id='stepsone'></a>
<ol>
<li>Elijamos un entero aleatorio $x$ entre $1$ y $N$ y calculamos el máximo común divisor $\text{mcd}(x,N)$ usando el algoritmo de Euclides.</li>
<li>Si $x$ y $N$ tienen algunos factores primos comunes, $\text{mcd}(x,N)$ será igual a $p$ o $q$. De lo contrario, $\text{gcd}(x,N) = 1$, lo que significa que $x$ y $N$ son coprimos. </li>
<li>Sea $r$ el período de $x \bmod N$ calculado por el algoritmo de búsqueda de períodos. Ahora hay que repetir los pasos anteriores con diferentes opciones aleatorias de $x$ hasta que $r$ sea par.</li>
<li>Ahora $p$ y $q$ se pueden encontrar calculando $\text{gcd}(x^{r/2} \pm 1, N)$ siempre que $x^{r/2} \neq \pm 1$.</li>
</ol>
</div>

Como ejemplo, consideremos $N = 15$. Veamos todos los valores de $1 < x < 15$ donde $x$ es coprimo con $15$:

|  $x$  |         $x^a \bmod 15$       | Period $r$ |$\text{gcd}(x^{r/2}-1,15)$|$\text{gcd}(x^{r/2}+1,15)$ | 
|:-----:|:----------------------------:|:----------:|:------------------------:|:-------------------------:|
|   2   | 1,2,4,8,1,2,4,8,1,2,4...     |      4     |             3            |             5             |
|   4   | 1,4,1,4,1,4,1,4,1,4,1...     |      2     |             3            |             5             |
|   7   | 1,7,4,13,1,7,4,13,1,7,4...   |      4     |             3            |             5             |
|   8   | 1,8,4,2,1,8,4,2,1,8,4...     |      4     |             3            |             5             |
|   11  | 1,11,1,11,1,11,1,11,1,11,1...|      2     |             5            |             3             |
|   13  | 1,13,4,7,1,13,4,7,1,13,4,... |      4     |             3            |             5             |
|   14  | 1,14,1,14,1,14,1,14,1,14,1,,,|      2     |             1            |             15            |

Como puede verse, cualquier valor de $x$ excepto $14$ devolverá los factores de $15$, es decir, $3$ y $5$. $14$ es un ejemplo del caso especial donde $(x^{r/2} + 1)$ o $(x^{r/2} - 1)$ es un múltiplo de $N$ y por lo tanto tenemos que probar con otro $x$.

En general, se puede demostrar que este caso especial ocurre con poca frecuencia, por lo que, en promedio, solo dos llamadas a al algoritmo de determinación del período son suficientes para factorizar $N$.

### Búsqueda cuántica del período  <a id='quantumperiodfinding'></a>

Primero describamos el algoritmo de búsqueda cuántica del período y luego veamos algunos de los pasos en detalle, antes de pasar a un ejemplo. Este algoritmo toma dos enteros coprimos, $x$ y $N$, y genera $r$, el período de $\mathcal{F}(a) = x^a\bmod N$.

<div class="alert alert-block alert-info"><a id='stepstwo'></a>
<ol>
<li> Elegir $T = 2^t$ tal que $N^2 \leq T \le 2N^2$. Inicializar dos registros de qubits, primero un registro de argumento con $t$ qubits y segundo un registro de función con $n = log_2 N$ qubits. Estos registros comienzan en el estado inicial:
$$\vert\psi_0\rangle = \vert 0 \rangle \vert 0 \rangle$$ </li>
<li> Aplicar una puerta de Hadamard en cada uno de los qubits en el registro de argumento para obtener una superposición de todos los números enteros de $0$ a $T$ con el mismo peso:
$$\vert\psi_1\rangle = \frac{1}{\sqrt{T}}\sum_{a=0}^{T-1}\vert a \rangle \vert 0 \rangle$$ </li>
<li> Implementar la función de exponenciación modular $x^a \bmod N$ en el registro de función, dando el estado:
$$\vert\psi_2\rangle = \frac{1}{\sqrt{T}}\sum_{a=0}^{T-1}\vert a \rangle \vert x^a \bmod N \rangle$$
Este $\vert\psi_2\rangle$ está altamente entrelazado y muestra paralelismo cuántico, es decir, la función entrelazó en paralelo todos los valores de entrada de 0 a $T$ con los valores correspondientes de $x^a \bmod N$, aunque la función se ejecutó una vez. </li>
<li> Realizar la transformada cuántica de Fourier en el registro de argumentos, lo que da como resultado el estado:
$$\vert\psi_3\rangle = \frac{1}{T}\sum_{a=0}^{T-1}\sum_{z=0}^{T-1}e^{(2\pi i)(az/T)}\vert z \rangle \vert x^a \bmod N \rangle$$
donde debido a la interferencia, solo los términos $\vert z \rangle$ con
$$z = qT/r $$
tienen una amplitud significativa donde $q$ es un número entero aleatorio que va de $0$ a $r-1$ y $r$ es el período de $\mathcal{F}(a) = x^a\bmod N$. </li>
<li> Medir el registro de argumento para obtener el resultado clásico $z$. Con una probabilidad razonable, la aproximación de fracción continua de $T / z$ será un múltiplo entero del período $r$. El algoritmo de Euclides se puede usar para encontrar $r$.</li>
</ol>
</div>

Se puede observar cómo se han utilizado el paralelismo cuántico y la interferencia constructiva para detectar y medir la periodicidad de la función de exponenciación modular. El hecho de que la interferencia facilite la medición de la periodicidad no debería ser una gran sorpresa. Ya que, los físicos utilizan habitualmente la dispersión de ondas electromagnéticas y las mediciones de interferencia para determinar la periodicidad de los objetos físicos, como las redes cristalinas. Asimismo, el algoritmo de Shor explota la interferencia para medir la periodicidad de los objetos aritméticos.

#### Exponenciación modular

La exponenciación modular, paso 3 anterior, que es la evaluación de $x^a \bmod N$ para $2^t$ valores de $a$ en paralelo, es la parte más exigente del algoritmo (la que genera los cuello de botella en la implementación real del circuito cuántico). Esta se puede realizar usando la siguiente identidad para la representación binaria de cualquier número entero: $x = x_{t-1}2^{t-1} + ... x_12^1+x_02^0$, donde $x_t$ son los dígitos binarios de $x$. De esto se sigue que:

$$ x^a \bmod N = x^{2^{(t-1)}a_{t-1}} ... x^{2a_1}x^{a_0} \bmod N = x^{2^{(t-1)}a_{t-1}}  ... [x^{2a_1}[x^{2a_0} \bmod N] \bmod N]  ... \bmod N $$

Esto significa que 1 se multiplica primero por $x^1 \bmod N$ si y solo si $a_0 = 1$, luego el resultado se multiplica por $x^2 \bmod N$ si y solo si $a_1 = 1$ y así sucesivamente, hasta que finalmente el resultado se multiplica por $x^{2^{(s-1)}}\bmod N$ si y solo si $a_{t-1} = 1$.

Por tanto, la exponenciación modular consta de $t$ multiplicaciones en serie módulo $N$, cada una de ellas controlada por el qubit $a_t$. Los valores $x,x^2,...,x^{2^{(t-1)}} \bmod N$ se pueden encontrar eficientemente en una computadora clásica mediante el algoritmo de cuadrados repetidos.

#### Transformada cuántica de Fourier

La transformada de Fourier ocurre en muchas versiones diferentes a lo largo de la computación clásica, en áreas que van desde el procesamiento de señales hasta la compresión de datos y la teoría de la complejidad. La transformada cuántica de Fourier (QFT), el paso 4 anterior, es la implementación cuántica de la transformada discreta de Fourier sobre las amplitudes de una función de onda.

La transformada de Fourier discreta clásica actúa sobre un vector $(x_0, ..., x_{N-1})$ y lo asigna al vector $(y_0, ..., y_{N-1})$ de acuerdo con la fórmula
$$y_k = \frac{1}{\sqrt{N}}\sum_{j=0}^{N-1}x_j\omega_N^{jk}$$
donde $\omega_N^{jk} = e^{2\pi i \frac{jk}{N}}$.

De manera similar, la transformada cuántica de Fourier actúa sobre un estado cuántico $\sum_{i=0}^{N-1} x_i \vert i \rangle$ y lo asigna al estado cuántico $\sum_{i=0}^{N -1} y_i \vert i \rangle$ según la fórmula
$$y_k = \frac{1}{\sqrt{N}}\sum_{j=0}^{N-1}x_j\omega_N^{jk}$$
con $\omega_N^{jk}$ definido como arriba. Hay que tener en cuenta que solo las amplitudes del estado se vieron afectadas por esta transformación.

Esto también se puede expresar como la aplicación:
$$\vert x \rangle \mapsto \frac{1}{\sqrt{N}}\sum_{y=0}^{N-1}\omega_N^{xy} \vert y \rangle$$

O la matriz unitaria:
$$ U_{QFT} = \frac{1}{\sqrt{N}} \sum_{x=0}^{N-1} \sum_{y=0}^{N-1} \omega_N^{xy } \vert y \rangle \langle x \vert$$

Como ejemplo, ya hemos visto la transformada cuántica de Fourier cuando $N = 2$, es el operador de Hadamard ($H$):
$$H = \frac{1}{\sqrt{2}}\begin{bmatrix} 1 & 1 \\ 1 & -1 \end{bmatrix}$$
Supongamos que tenemos el estado de qubit único $\alpha \vert 0 \rangle + \beta \vert 1 \rangle$, si aplicamos el operador $H$ a este estado, obtenemos el nuevo estado:
$$\frac{1}{\sqrt{2}}(\alpha + \beta) \vert 0 \rangle + \frac{1}{\sqrt{2}}(\alpha - \beta) \vert 1 \rangle \equiv \tilde{\alpha}\vert 0 \rangle + \tilde{\beta}\vert 1 \rangle$$
Se observa cómo la puerta de Hadamard realiza la transformada discreta de Fourier para $N = 2$ en las amplitudes del estado.

Entonces, ¿cómo será la transformada cuántica de Fourier para N más grandes? Construyamos un circuito para $N=2^n$, $QFT_N$ actuando en el estado $\vert x \rangle = \vert x_1...x_n \rangle$ donde $x_1$ es el bit más significativo.

$$QFT_N\vert x \rangle = \frac{1}{\sqrt{N}} \sum_{y=0}^{N-1}\omega_N^{xy} \vert y \rangle = \frac{1}{\sqrt{N}} \sum_{y=0}^{N-1} e^{2 \pi i xy / 2^n} \vert y \rangle$$ ya que $\omega_N^{xy} = e^{2\pi i \frac{xy}{N}}$ y $N = 2^n$
$$= \frac{1}{\sqrt{N}} \sum_{y=0}^{N-1} e^{2 \pi i \left(\sum_{k=1}^n y_k/2^k\right) x} \vert y_1 ... y_n \rangle$$ lo reescribimos en notación fraccional binaria $y = y_1...y_k, y/2^n = \sum_{k=1}^n y_k/2^k$
$$= \frac{1}{\sqrt{N}} \sum_{y=0}^{N-1} \prod_{k=0}^n e^{2 \pi i x y_k/2^k } \vert y_1 ... y_n \rangle$$ después de expandir la exponencial de una suma en un producto de exponenciales
$$= \frac{1}{\sqrt{N}} \bigotimes_{k=1}^n  \left(\vert0\rangle + e^{2 \pi i x /2^k } \vert1\rangle \right)$$ después de reordenar la suma y los productos, y expandir} $$= \frac{1}{\sqrt{N}} \left(\vert0\rangle + e^{2 \pi i[0.x_n]} \vert1\rangle\right) \otimes...\otimes  \left(\vert0\rangle + e^{2 \pi i[0.x_1.x_2...x_{n-1}.x_n]} \vert1\rangle\right)$$ como $e^{2 \pi i x/2^k} = e^{2 \pi i[0.x_k...x_n]}$

Esta es una forma muy útil de QFT para $N=2^n$ ya que solo el último qubit depende del valor de todos los demás qubits de entrada y cada bit adicional depende cada vez menos de los qubits de entrada. Además, hay que tener en cuenta que $e^{2 \pi i.0.x_n}$ es $+1$ o $-1$, se asemeja a la transformación de Hadamard.

Antes de crear el código de circuito para el general $N=2^n$, veamos $N=8,n=3$:
$$QFT_8\vert x_1x_2x_3\rangle = \frac{1}{\sqrt{8}} \left(\vert0\rangle + e^{2 \pi i[0.x_3]} \vert1\rangle\right) \otimes \left(\vert0\rangle + e^{2 \pi i[0.x_2.x_3]} \vert1\rangle\right) \otimes \left(\vert0\rangle + e^{2 \pi i[0 .x_1.x_2.x_3]} \vert1\rangle\right) $$

Los pasos para crear el circuito para $\vert y_1y_2x_3\rangle = QFT_8\vert x_1x_2x_3\rangle$, recordamos que la puerta de rotación de fase controlada $CU_1$, sería:
1. Aplicar una puerta Hadamard a $\vert x_3 \rangle$, dando el estado $\frac{1}{\sqrt{2}}\left(\vert0\rangle + e^{2 \pi i.0.x_3} \vert1\rangle\right) = \frac{1}{\sqrt{2}}\left(\vert0\rangle + (-1)^{x_3} \vert1\rangle\right)$
2. Aplicar una puerta Hadamard a $\vert x_2 \rangle$, luego dependiendo de $k_3$ (antes de la puerta de Hadamard) un $CU_1(\frac{\pi}{2})$, dando el estado $\frac{1 }{\sqrt{2}}\left(\vert0\rangle + e^{2 \pi i[0.x_2.x_3]} \vert1\rangle\right)$.
3. Apliar una puerta Hadamard a $\vert x_1 \rangle$, luego $CU_1(\frac{\pi}{2})$ dependiendo de $k_2$, y $CU_1(\frac{\pi}{4})$ dependiendo de $k_3$.
4. Medir los bits en orden inverso, es decir, $y_3 = x_1, y_2 = x_2, y_1 = y_3$.

#### Ejemplo

Factoricemos $N = 21$ con el coprimo $x=2$, siguiendo los pasos anteriores del algoritmo de búsqueda del período cuántico, nos debería devolver $r = 6$. 

1. Elegimos $T = 2^t$ tal que $N^2 \leq T \le 2N^2$. Para $N = 21$, el valor más pequeño de $t$ es 9, lo que significa que $T = 2^t = 512$. Inicializamos dos registros de qubits, primero un registro de argumento con $t = 9$ qubits y segundo un registro de función con $n = log_2 N = 5$ qubits:
$$\vert\psi_0\rangle = \vert 0 \rangle \vert 0 \rangle$$

2. Aplicamos una puerta de Hadamard a cada uno de los qubits en el registro de argumento:
$$\vert\psi_1\rangle = \frac{1}{\sqrt{T}}\sum_{a=0}^{T-1}\vert a \rangle \vert 0 \rangle = \frac{1} {\sqrt{512}}\sum_{a=0}^{511}\vert a \rangle \vert 0 \rangle$$

3. Implementamos la función de exponenciación modular $x^a \bmod N$ en el registro de función:

$$\vert\psi_2\rangle  =  \frac{1}{\sqrt{T}}\sum_{a=0}^{T-1}\vert a \rangle \vert x^a \bmod N \rangle = \frac{1}{\sqrt{512}}\sum_{a=0}^{511}\vert a \rangle \vert 2^a \bmod 21 \rangle$$
$$= \frac{1}{\sqrt{512}} \bigg( \vert 0 \rangle \vert 1 \rangle + \vert 1 \rangle \vert 2 \rangle + \vert 2 \rangle \vert 4 \rangle + \vert 3 \rangle \vert 8 \rangle + \vert 4 \rangle \vert 16 \rangle + \vert 5 \rangle \vert 11 \rangle + \vert 6 \rangle \vert 1 \rangle + \vert 7 \rangle \vert 2 \rangle + \vert 8 \rangle \vert 4 \rangle + \vert 9 \rangle \vert 8 \rangle + \vert 10 \rangle \vert 16 \rangle + \vert 11 \rangle \vert 11 \rangle \vert 12 \rangle \vert 1 \rangle + \ldots \bigg)$$
Observamos que la expresión anterior tiene el siguiente patrón: los estados del segundo registro de cada "columna" son los mismos. Por tanto podemos reordenar los términos para recoger el segundo registro:

$$\vert\psi_2\rangle = \frac{1}{\sqrt{512}} \bigg[ \big(\vert 0 \rangle + \vert 6 \rangle + \vert 12 \rangle \ldots + \vert 504 \rangle + \vert 510 \rangle \big) \vert 1 \rangle + \big(\vert 1 \rangle + \vert 7 \rangle + \vert 13 \rangle \ldots + \vert 505 \rangle + \vert 511 \rangle \big) \vert 2 \rangle + \big(\vert 2 \rangle + \vert 8 \rangle + \vert 14 \rangle \ldots + \vert 506 \rangle +  \big) \vert 4 \rangle$$

$$+ \big(\vert 3 \rangle + \vert 9 \rangle + \vert 15 \rangle \ldots + \vert 507 \rangle +  \big) \vert 8 \rangle + \big(\vert 4 \rangle + \vert 10 \rangle + \vert 16 \rangle \ldots + \vert 508 \rangle + \big) \vert 16 \rangle + \big(\vert 5 \rangle + \vert 11 \rangle + \vert 17 \rangle \ldots + \vert 509 \rangle +  \big) \vert 11 \rangle \bigg]$$

4. Para simplificar las siguientes ecuaciones, mediremos el registro de función antes de realizar una transformada cuántica de Fourier en el registro del argumento. Esto producirá uno de los siguientes números con la misma probabilidad: $\{1,2,4,6,8,16,11\}$. Supongamos que el resultado de la medición fue $2$, entonces:
$$\vert\psi_3\rangle = \frac{1}{\sqrt{86}}(\vert 1 \rangle + \vert 7 \rangle + \vert 13 \rangle \ldots + \vert 505 \rangle + \vert 511 \rangle) \vert 2 \rangle $$
No importa cuál sea el resultado de la medición; lo que importa es el patrón periódico. El periodo de los estados del primer registro es la solución al problema y la transformada cuántica de Fourier puede revelar el valor del periodo.

5. Realizamos una transformada cuántica de Fourier en el registro de argumento:
$$\vert\psi_4\rangle = QFT(\vert\psi_3\rangle) = QFT(\frac{1}{\sqrt{86}}\sum_{a=0}^{85}\vert 6a+1 \rangle)\vert 2 \rangle = \frac{1}{\sqrt{512}}\sum_{j=0}^{511}\bigg(\big[ \frac{1}{\sqrt{86}}\sum_{a=0}^ {85} e^{-2 \pi i \frac{6ja}{512}} \big] e^{-2\pi i\frac{j}{512}}\vert j \rangle \bigg)\vert 2 \rangle$$

6. Medimos el registro de argumento. La probabilidad de medir un resultado $j$ es:
$$ \rm{Probabilidad}(j) = \frac{1}{512 \times 86} \bigg\vert \sum_{a=0}^{85}e^{-2 \pi i \frac{6ja} {512}} \bigg\vert^2$$
Esto alcanza un máximo en $j=0,85,171,256,341,427$. Supongamos que el resultado de la medición nos dió $j = 85$, ahora usando la aproximación de fracción continua de $\frac{512}{85}$, obtenemos $r=6$, como se esperaba.
