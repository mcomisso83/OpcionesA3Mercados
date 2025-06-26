/**
 * Calcula el valor de una opción europea y sus griegas con el modelo de Black-Scholes Generalizado.
 * @param {string} parametro Parámetro a calcular: prima, delta, gamma, gammap, theta, vega, rho.
 * @param {number} s Precio spot del subyacente.
 * @param {number} k Precio de ejercicio.
 * @param {number} T Tiempo hasta el vencimiento (en años).
 * @param {number} r Tasa libre de riesgo.
 * @param {number} q Tasa de dividendos o costo de carry.
 * @param {number} sigma Volatilidad del subyacente.
 * @param {string} CallPut "call" o "put".
 * @param {string} tipo_subyacente "accion" o "futuro".
 * @param {string} estilo_liquidacion "futures_style" o "equity_style". Solo aplica si tipo_subyacente = "futuro".
 * @return El valor del parámetro solicitado.
 * @customfunction
 */
function BlackScholes(parametro, s, k, T, r, q, sigma, CallPut, tipo_subyacente, estilo_liquidacion) {
  parametro = parametro.toLowerCase();
  CallPut = CallPut.toLowerCase();
  tipo_subyacente = tipo_subyacente.toLowerCase();
  estilo_liquidacion = estilo_liquidacion ? estilo_liquidacion.toLowerCase() : "";

  const normaldf = x => Math.exp(-(x ** 2) / 2) / Math.sqrt(2 * Math.PI);
  const N = x => 0.5 * (1 + erf(x / Math.SQRT2));

  let b;

  if (tipo_subyacente === "futuro") {
    if (estilo_liquidacion !== "futures_style" && estilo_liquidacion !== "equity_style") {
      throw new Error("Indicar estilo de liquidación: 'futures_style' o 'equity_style'");
    }

    b = 0;
    if (estilo_liquidacion === "futures_style") {
      r = 0;
    }
  } else if (tipo_subyacente === "accion") {
    b = r - q;
  } else {
    throw new Error("tipo_subyacente debe ser 'accion' o 'futuro'.");
  }

  if (T <= 0) {
    const payoff = CallPut === "call" ? Math.max(s - k, 0) : Math.max(k - s, 0);
    const delta = CallPut === "call" ? (s > k ? 1 : 0) : (s < k ? -1 : 0);
    return parametro === "prima" ? payoff : parametro === "delta" ? delta : 0;
  }

  const d1 = (Math.log(s / k) + (b + sigma ** 2 / 2) * T) / (sigma * Math.sqrt(T));
  const d2 = d1 - sigma * Math.sqrt(T);
  const df_b = Math.exp((b - r) * T);
  const df_r = Math.exp(-r * T);

  let value = {};

  if (CallPut === "call") {
    value.prima = s * df_b * N(d1) - k * df_r * N(d2);
    value.delta = df_b * N(d1);
    value.theta = ((-s * df_b * normaldf(d1) * sigma) / (2 * Math.sqrt(T))
                 + (b - r) * s * df_b * N(d1)
                 - r * k * df_r * N(d2));
    value.gamma = (df_b * normaldf(d1)) / (s * sigma * Math.sqrt(T));
    value.gammap = value.gamma * s / 100;
    value.vega = s * df_b * normaldf(d1) * Math.sqrt(T) / 100;
    value.rho = tipo_subyacente === "futuro"
                ? (estilo_liquidacion === "equity_style" ? -T * value.prima / 100 : 0)
                : T * k * df_r * N(d2) / 100;
  } else {
    value.prima = k * df_r * N(-d2) - s * df_b * N(-d1);
    value.delta = -df_b * N(-d1);
    value.theta = ((-s * df_b * normaldf(d1) * sigma) / (2 * Math.sqrt(T))
                 + (b - r) * s * df_b * N(-d1)
                 + r * k * df_r * N(-d2));
    value.gamma = (df_b * normaldf(d1)) / (s * sigma * Math.sqrt(T));
    value.gammap = value.gamma * s / 100;
    value.vega = s * df_b * normaldf(d1) * Math.sqrt(T) / 100;
    value.rho = tipo_subyacente === "futuro"
                ? (estilo_liquidacion === "equity_style" ? -T * value.prima / 100 : 0)
                : -T * k * df_r * N(-d2) / 100;
  }

  return value[parametro] !== undefined ? value[parametro] : "Parámetro inválido";
}

/**
 * Estima la volatilidad implícita bajo Black-Scholes mediante bisección.
 * @param {number} s Precio del subyacente.
 * @param {number} k Strike.
 * @param {number} T Tiempo hasta vencimiento.
 * @param {number} r Tasa libre de riesgo.
 * @param {number} q Tasa de dividendos.
 * @param {number} prima Prima observada.
 * @param {string} CallPut "call" o "put"
 * @param {string} tipo_subyacente "accion" o "futuro"
 * @param {string} estilo_liquidacion "equity_style" o "futures_style". No aplica si tipo_subyacente = "accion".
 * @return Volatilidad implícita estimada.
 * @customfunction
 */
function BlackScholesVolImpl(s, k, T, r, q, prima, CallPut, tipo_subyacente, estilo_liquidacion) {
  CallPut = CallPut.toLowerCase();
  tipo_subyacente = tipo_subyacente.toLowerCase();
  estilo_liquidacion = estilo_liquidacion ? estilo_liquidacion.toLowerCase() : "";

  if (CallPut === "futuro") return "NA";
  if (prima <= 0 || T <= 0) return 0;

  let high = 4.0, low = 0.0, epsilon = 0.0001, counter = 0;
  while ((high - low) > epsilon) {
    counter++;
    if (counter > 100) return "NA";
    let mid = (high + low) / 2;
    let val = BlackScholes("prima", s, k, T, r, q, mid, CallPut, tipo_subyacente, estilo_liquidacion);
    if (val > prima) {
      high = mid;
    } else {
      low = mid;
    }
  }
  return (high + low) / 2;
}

// Función de error (para distribución normal acumulada)
function erf(x) {
  let a1 =  0.254829592, a2 = -0.284496736, a3 =  1.421413741;
  let a4 = -1.453152027, a5 =  1.061405429, p  =  0.3275911;
  let sign = x < 0 ? -1 : 1;
  x = Math.abs(x);
  let t = 1.0 / (1.0 + p * x);
  let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-(x * x));
  return sign * y;
}
