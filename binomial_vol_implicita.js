/**
 * Estima la volatilidad implícita bajo el modelo binomial.
 * @param {number} s Precio actual del subyacente.
 * @param {number} k Precio de ejercicio.
 * @param {number} T Tiempo hasta vencimiento (en años).
 * @param {number} r Tasa libre de riesgo.
 * @param {number} q Tasa de dividendos o acarreo.
 * @param {number} prima Prima observada.
 * @param {string} CallPut Tipo de opción: "call" o "put".
 * @param {string} EuroAmer Estilo de opción: "europea" o "americana".
 * @param {number} n Número de pasos del árbol binomial.
 * @param {string} tipo_subyacente Tipo de subyacente: "accion" o "futuro".
 * @param {string} estilo_liquidacion "futures_style" o "equity_style". Solo aplica si tipo_subyacente = "futuro".
 * @return Volatilidad implícita estimada o "NA" si no converge.
 * @customfunction
 */
function BinomialVolImpl(s, k, T, r, q, prima, CallPut, EuroAmer, n, tipo_subyacente, estilo_liquidacion) {
  CallPut = CallPut.toLowerCase();
  EuroAmer = EuroAmer.toLowerCase();
  tipo_subyacente = tipo_subyacente.toLowerCase();
  estilo_liquidacion = estilo_liquidacion ? estilo_liquidacion.toLowerCase() : "";

  if (CallPut === "futuro") return "NA";
  if (prima <= 0 || T <= 0 || n <= 0) return 0;

  if (tipo_subyacente === "futuro") {
    if (!["futures_style", "equity_style", "matba_rofex_style"].includes(estilo_liquidacion)) {
      throw new Error("Indicar estilo de liquidación: 'futures_style', 'equity_style' o 'matba_rofex_style'");
    }
  }

  let high = 4.0;
  let low = 0.0;
  const epsilon = 0.0001;
  let counter = 0;

  while ((high - low) > epsilon) {
    counter++;
    if (counter > 100) return "NA";

    let midSigma = (high + low) / 2;
    let val = BinomialOpVal("prima", s, k, T, r, q, midSigma, CallPut, EuroAmer, n, tipo_subyacente, estilo_liquidacion);

    if (val > prima) {
      high = midSigma;
    } else {
      low = midSigma;
    }
  }

  return (high + low) / 2;
}
