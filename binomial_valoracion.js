/**
 * Calcula el valor de una opción y sus griegas utilizando el modelo binomial.
 * @param {string} parametro Parámetro a calcular: prima, delta, gamma, theta, vega, rho.
 * @param {number} s Precio del subyacente.
 * @param {number} k Precio de ejercicio.
 * @param {number} T Tiempo hasta vencimiento (en años).
 * @param {number} r Tasa libre de riesgo.
 * @param {number} q Tasa de dividendos o acarreo.
 * @param {number} sigma Volatilidad.
 * @param {string} CallPut Tipo de opción: "call" o "put".
 * @param {string} EuroAmer Estilo de opción: "europea" o "americana".
 * @param {number} n Número de pasos del árbol binomial.
 * @param {string} tipo_subyacente Tipo de subyacente: "accion" o "futuro".
 * @param {string} estilo_liquidacion "futures_style", "equity_style" o "matba_rofex_style". Solo aplica si tipo_subyacente = "futuro".
 * @return Valor del parámetro solicitado.
 * @customfunction
 */
function BinomialOpVal(parametro, s, k, T, r, q, sigma, CallPut, EuroAmer, n, tipo_subyacente, estilo_liquidacion) {
  if (!Number.isInteger(n) || n <= 1) throw new Error("El número de pasos 'n' debe ser un entero mayor que 1.");

  CallPut = CallPut.toLowerCase();
  EuroAmer = EuroAmer.toLowerCase();
  tipo_subyacente = tipo_subyacente.toLowerCase();
  estilo_liquidacion = estilo_liquidacion ? estilo_liquidacion.toLowerCase() : "";

  const z = CallPut === "call" ? 1 : -1;
  let b;

  if (tipo_subyacente === "futuro") {
    if (!["futures_style", "equity_style", "matba_rofex_style"].includes(estilo_liquidacion)) {
      throw new Error("Indicar estilo de liquidación: 'futures_style', 'equity_style' o 'matba_rofex_style'");
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
    let payoff = Math.max(0, z * (s - k));
    if (parametro === "prima") return payoff;
    if (parametro === "delta") return z * (s !== k ? 1 : 0);
    return 0;
  }

  const dt = T / n;
  const u = Math.exp(sigma * Math.sqrt(dt));
  const d = 1 / u;
  const p = (Math.exp(b * dt) - d) / (u - d);
  const Df = Math.exp(-r * dt);
  const interestFactor = Math.exp(r * dt) - 1;

  let OptionValue = Array(n + 1).fill(0);
  let Delta = 0, Gamma = 0, Theta = 0;
  let C00 = 0, C21 = 0;

  for (let i = 0; i <= n; i++) {
    let ST = s * Math.pow(u, i) * Math.pow(d, n - i);
    OptionValue[i] = Math.max(0, z * (ST - k));
  }

  for (let j = n - 1; j >= 0; j--) {
    for (let i = 0; i <= j; i++) {
      let ST = s * Math.pow(u, i) * Math.pow(d, j - i);
      let intrinsic = Math.max(0, z * (ST - k));
      let futureVal = p * OptionValue[i + 1] + (1 - p) * OptionValue[i];

      let expected;
      if (tipo_subyacente === "futuro" && estilo_liquidacion === "matba_rofex_style") {
        expected = Df * (intrinsic * interestFactor + futureVal);
      } else {
        expected = Df * futureVal;
      }

      OptionValue[i] = EuroAmer === "americana"
        ? Math.max(intrinsic, expected)
        : expected;
    }

    if (j === 2) {
      let Su = s * Math.pow(u, 2);
      let Sm = s * u * d;
      let Sd = s * Math.pow(d, 2);
      let Cuu = OptionValue[2];
      let Cu = OptionValue[1];
      let Cd = OptionValue[0];
      Gamma = ((Cuu - Cu) / (Su - Sm) - (Cu - Cd) / (Sm - Sd)) / ((Su - Sd) / 2);
      C21 = OptionValue[1];
    }

    if (j === 1) {
      let Su = s * u;
      let Sd = s * d;
      Delta = (OptionValue[1] - OptionValue[0]) / (Su - Sd);
    }

    if (j === 0) {
      C00 = OptionValue[0];
    }
  }

  Theta = (C21 - C00) / (2 * dt);
  const Prima = OptionValue[0];

  if (parametro === "vega") {
    const change = 0.01;
    const up = BinomialOpVal("prima", s, k, T, r, q, sigma + change, CallPut, EuroAmer, n, tipo_subyacente, estilo_liquidacion);
    const down = BinomialOpVal("prima", s, k, T, r, q, sigma - change, CallPut, EuroAmer, n, tipo_subyacente, estilo_liquidacion);
    return (up - down) / 2;
  }

  if (parametro === "rho") {
    const change = 0.01;
    const up = BinomialOpVal("prima", s, k, T, r + change, q, sigma, CallPut, EuroAmer, n, tipo_subyacente, estilo_liquidacion);
    const down = BinomialOpVal("prima", s, k, T, r - change, q, sigma, CallPut, EuroAmer, n, tipo_subyacente, estilo_liquidacion);
    return (up - down) / 2;
  }

  switch (parametro) {
    case "prima": return Prima;
    case "delta": return Delta;
    case "gamma": return Gamma;
    case "theta": return Theta;
    default: return "Parámetro no reconocido.";
  }
}
