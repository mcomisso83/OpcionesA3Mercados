# OpcionesA3Mercados
Scripts de valoración de opciones sobre futuros con estilos de liquidación Matba Rofex, equity y futures-style, implementados en Google Apps Script.
# OpcionesMatbaRofex

Este repositorio contiene scripts desarrollados en Google Apps Script para la valoración de opciones sobre futuros en diferentes estilos de liquidación:

- **Estilo Matba Rofex**
- **Equity-style**
- **Futures-style**

Los modelos incluyen:
- Árbol binomial adaptado (con griegas)
- Estimación de volatilidad implícita
- Modelo generalizado de Black-Scholes

## Archivos

- `binomial_valoracion.js`: Valoración binomial con soporte para estilos de liquidación.
- `binomial_volatilidad_implicita.js`: Cálculo de volatilidad implícita por bisección.
- `black_scholes_generalizado.js`: Implementación del modelo de Black-Scholes con griegas.

> Estos scripts están pensados para integrarse con Google Sheets mediante Apps Script.

## Licencia

MIT License
