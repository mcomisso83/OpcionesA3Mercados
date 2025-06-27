# OpcionesA3Mercados

Scripts de valoración de opciones sobre futuros con estilos de liquidación **Matba Rofex**, **equity-style** y **futures-style**, implementados en Google Apps Script.

## Estilos de liquidación soportados

- **Matba Rofex-style**: ajuste diario sobre el valor intrínseco.
- **Equity-style**: prima pagada al inicio, sin MtM diario.
- **Futures-style**: liquidación al vencimiento con MtM diario sobre el valor de la opción.

## Modelos incluidos

- ✅ Árbol binomial recombinante con griegas (delta, gamma, theta, vega, rho)
- ✅ Estimación de volatilidad implícita (binomial y Black-Scholes)
- ✅ Modelo generalizado de Black-Scholes (adaptado a estilos sobre futuros)

---

## Archivos disponibles

### 🔹 Implementación modular (recomendada para desarrollo)

- `black_scholes_generalizado.js`: Modelo de Black-Scholes con griegas.
- `binomial_valoracion.js`: Valoración binomial con soporte para los tres estilos.
- `binomial_volatilidad_implicita.js`: Cálculo de volatilidad implícita bajo el modelo binomial.

Estos archivos pueden copiarse individualmente dentro del editor de Google Apps Script.

---

### 🔸 Implementación todo-en-uno (recomendada para usuarios finales)

- `OpcionesA3Mercados_TodoEnUno.gs` (o `.txt`): contiene todas las funciones unificadas en un solo script.

Ideal para pegar directamente en un solo proyecto de Google Sheets y comenzar a usar las funciones `=BlackScholes(...)`, `=BinomialOpVal(...)`, etc.

---

## Cómo usar en Google Sheets

1. Abrí tu hoja de cálculo.
2. Ir a **Extensiones > Apps Script**.
3. Pegá el contenido de uno o varios archivos `.js` o el archivo unificado.
4. Guardá y autorizá los permisos.
5. ¡Listo! Ya podés usar las funciones personalizadas en tus celdas.

---

## Licencia

MIT License

