// === Black-Scholes ===

/**
 * Calcula el valor de una opción europea y sus griegas con el modelo de Black-Scholes Generalizado.
 * @param {string} parametro Parámetro a calcular: "prima", "delta", "gamma", "gammap", "theta", "vega", "rho".
 * @param {number} s Precio spot del subyacente.
 * @param {number} k Precio de ejercicio.
 * @param {number} T Tiempo hasta el vencimiento (en años).
 * @param {number} r Tasa libre de riesgo.
 * @param {number} q Tasa de dividendos o costo de carry.
 * @param {number} sigma Volatilidad del subyacente (valor positivo, en forma decimal).
 * @param {string} CallPut "call" o "put".
 * @param {string} tipo_subyacente "accion" o "futuro".
 * @param {string} estilo_liquidacion "futures_style" o "equity_style". Solo aplica si tipo_subyacente = "futuro".
 * @return El valor del parámetro solicitado (número). Si el parámetro es inválido, devuelve un mensaje de error.
 * @customfunction
 */
function BlackScholes(parametro, s, k, T, r, q, sigma, CallPut, tipo_subyacente, estilo_liquidacion) {
  parametro = parametro.toLowerCase();
  CallPut = CallPut.toLowerCase();
  tipo_subyacente = tipo_subyacente.toLowerCase();
  estilo_liquidacion = estilo_liquidacion ? estilo_liquidacion.toLowerCase() : "";

  const normaldf = x => Math.exp(-(x ** 2) / 2) / Math.sqrt(2 * Math.PI);

  const N = function(x) {
    const a1 = 0.319381530, a2 = -0.356563782, a3 = 1.781477937,
          a4 = -1.821255978, a5 = 1.330274429;
    const k = 1.0 / (1.0 + 0.2316419 * Math.abs(x));
    const poly = ((((a5 * k + a4) * k + a3) * k + a2) * k + a1) * k;
    const approx = 1.0 - normaldf(x) * poly;
    return x >= 0 ? approx : 1 - approx;
  };

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
 * @param {number} T Tiempo hasta vencimiento (en años).
 * @param {number} r Tasa libre de riesgo.
 * @param {number} q Tasa de dividendos o costo de carry.
 * @param {number} prima Prima observada (valor de mercado).
 * @param {string} CallPut "call" o "put".
 * @param {string} tipo_subyacente "accion" o "futuro".
 * @param {string} estilo_liquidacion "equity_style" o "futures_style". No aplica si tipo_subyacente = "accion".
 * @return Volatilidad implícita estimada. Si no converge, devuelve "NA".
 * @customfunction
 */
function BlackScholesVolImpl(s, k, T, r, q, prima, CallPut, tipo_subyacente, estilo_liquidacion) {
  CallPut = CallPut.toLowerCase();
  tipo_subyacente = tipo_subyacente.toLowerCase();
  estilo_liquidacion = estilo_liquidacion ? estilo_liquidacion.toLowerCase() : "";

  if (tipo_subyacente !== "accion" && tipo_subyacente !== "futuro") return "NA";
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


// === BinomialOpVal ===

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

// === BinomialVolImpl ===

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
 * @param {string} estilo_liquidacion "futures_style", "equity_style" o "matba_rofex_style". Solo aplica si tipo_subyacente = "futuro".
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

// === Bjerksund & Stensland (1993) ===

/**
 * Calcula la prima o las griegas de una opción americana usando el modelo de Bjerksund y Stensland (1993).
 * @param {string} parametro "prima" | "delta" | "gamma" | "theta" | "vega" | "rho"
 * @param {number} s Precio spot del subyacente
 * @param {number} k Strike
 * @param {number} T Tiempo a vencimiento (años)
 * @param {number} r Tasa libre de riesgo (decimal)
 * @param {number} q Tasa de dividendos/carry (decimal)
 * @param {number} sigma Volatilidad (decimal)
 * @param {string} CallPut "call" | "put"
 * @param {string} tipo_subyacente "accion" | "futuro"
 * @return {number}
 * @customfunction
 */
function BS93Amer(parametro, s, k, T, r, q, sigma, CallPut, tipo_subyacente) {
  parametro = String(parametro).toLowerCase();
  CallPut = String(CallPut).toLowerCase();
  tipo_subyacente = String(tipo_subyacente).toLowerCase();

  // b inicial (solo se usa en prima/delta/gamma/vega/theta; para rho se recalcula con r±h)
  let b;
  if (tipo_subyacente === "accion") {
    b = r - q;
  } else if (tipo_subyacente === "futuro") {
    b = 0;
  } else {
    throw new Error("tipo_subyacente debe ser 'accion' o 'futuro'");
  }

  const h = 0.01; // paso para diferencias finitas (1%)

  function f(s_, k_, T_, r_, b_, sigma_) {
    return BS93AmericanOption(CallPut, s_, k_, T_, r_, b_, sigma_);
  }

  switch (parametro) {
    case "prima":
      return f(s, k, T, r, b, sigma);

    case "delta":
      return (f(s + h, k, T, r, b, sigma) - f(s - h, k, T, r, b, sigma)) / (2 * h);

    case "gamma":
      return (
        f(s + h, k, T, r, b, sigma) - 2 * f(s, k, T, r, b, sigma) + f(s - h, k, T, r, b, sigma)
      ) / (h * h);

    case "theta": {
      // centrado en T para mejor precisión (resultado por año)
      const hT = Math.min(h, Math.max(1e-6, 0.5 * T)); // evita T-h negativo
      return (f(s, k, T - hT, r, b, sigma) - f(s, k, T + hT, r, b, sigma)) / (2 * hT);
    }
    
    case "vega":
      // por 1% de sigma
      return (f(s, k, T, r, b, sigma + h) - f(s, k, T, r, b, sigma - h)) / (2 * h) / 100;

    case "rho": {
      // ∂V/∂r con q constante:
      // si es "accion", b = r - q debe moverse con r; si es "futuro", b=0 constante
      const priceWithR = (rLocal) => {
        const bLocal =
          tipo_subyacente === "accion"
            ? (rLocal - q)
            : 0; // "futuro"
        return BS93AmericanOption(CallPut, s, k, T, rLocal, bLocal, sigma);
      };
      // derivada centrada con SIGNO CORRECTO; resultado por 1% de tasa
      return (priceWithR(r + h) - priceWithR(r - h)) / (2 * h) / 100;
    }

    default:
      return "Parámetro inválido";
  }
}

function BS93AmericanOption(optionType, S, X, T, r, b, sigma) {
  if (optionType.toLowerCase() === "put") {
    // Transformación put→call en el marco BS93
    return BS93AmericanCall(X, S, T, r - b, -b, sigma);
  } else if (optionType.toLowerCase() === "call") {
    return BS93AmericanCall(S, X, T, r, b, sigma);
  } else {
    throw new Error("El tipo de opción debe ser 'call' o 'put'");
  }
}

function BS93AmericanCall(S, X, T, r, b, sigma) {
  if (b >= r) return BS_EuropeanCall(S, X, T, r, b, sigma);

  const sig2 = sigma * sigma;
  const beta = (0.5 - b / sig2) + Math.sqrt(Math.pow(b / sig2 - 0.5, 2) + 2 * r / sig2);
  const BInfinity = (beta / (beta - 1)) * X;
  const B0 = Math.max(X, (r / (r - b)) * X);
  const hT = -(b * T + 2 * sigma * Math.sqrt(T)) * (B0 / (BInfinity - B0));
  const I = B0 + (BInfinity - B0) * (1 - Math.exp(hT));
  const alpha = (I - X) * Math.pow(I, -beta);

  if (S >= I) return S - X;

  return (
    alpha * Math.pow(S, beta)
    - alpha * phi(S, T, beta, I, I, r, b, sigma)
    + phi(S, T, 1, I, I, r, b, sigma)
    - phi(S, T, 1, X, I, r, b, sigma)
    - X * phi(S, T, 0, I, I, r, b, sigma)
    + X * phi(S, T, 0, X, I, r, b, sigma)
  );
}

function phi(S, T, gamma, H, I, r, b, sigma) {
  const sig2 = sigma * sigma;
  const lambda = -r + gamma * b + 0.5 * gamma * (gamma - 1) * sig2;
  const d = -(Math.log(S / H) + (b + (gamma - 0.5) * sig2) * T) / (sigma * Math.sqrt(T));
  const kappa = (2 * b / sig2) + (2 * gamma - 1);
  return Math.exp(lambda * T) * Math.pow(S, gamma) *
    (normCDF(d) - Math.pow(I / S, kappa) * normCDF(d - 2 * Math.log(I / S) / (sigma * Math.sqrt(T))));
}

function BS_EuropeanCall(S, X, T, r, b, sigma) {
  const d1 = (Math.log(S / X) + (b + 0.5 * sigma * sigma) * T) / (sigma * Math.sqrt(T));
  const d2 = d1 - sigma * Math.sqrt(T);
  return S * Math.exp((b - r) * T) * normCDF(d1) - X * Math.exp(-r * T) * normCDF(d2);
}

function normCDF(x) {
  return (1 + erf(x / Math.sqrt(2))) / 2;
}

function erf(x) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x);
  const t = 1 / (1 + p * x);
  const y = 1 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * Math.exp(-x * x);
  return sign * y;
}


/**
 * Estima la volatilidad implícita bajo el modelo de Bjerksund & Stensland (1993) mediante el método de bisección.
 * @param {number} s Precio del subyacente.
 * @param {number} k Strike.
 * @param {number} T Tiempo hasta vencimiento (en años).
 * @param {number} r Tasa libre de riesgo.
 * @param {number} q Tasa de dividendos o costo de carry.
 * @param {number} prima Prima observada (valor de mercado).
 * @param {string} CallPut "call" o "put".
 * @param {string} tipo_subyacente "accion" o "futuro".
 * @return Volatilidad implícita estimada. Si no converge, devuelve "NA".
 * @customfunction
 */
function BS93VolImpl(s, k, T, r, q, prima, CallPut, tipo_subyacente) {
  CallPut = CallPut.toLowerCase();
  tipo_subyacente = tipo_subyacente.toLowerCase();

  if (tipo_subyacente !== "accion" && tipo_subyacente !== "futuro") return "NA";
  if (prima <= 0 || T <= 0) return 0;

  let low = 0.0001, high = 3.0, epsilon = 0.0001, maxIter = 100;
  let counter = 0;

  while ((high - low) > epsilon) {
    counter++;
    if (counter > maxIter) return "NA";
    let mid = (high + low) / 2;
    let val = BS93Amer("prima", s, k, T, r, q, mid, CallPut, tipo_subyacente);
    if (val > prima) {
      high = mid;
    } else {
      low = mid;
    }
  }

  return (high + low) / 2;
}
