/**
 * Calcula la prima o las griegas de una opción americana usando el modelo de Bjerksund y Stensland (1993).
 * @param {string} parametro Parámetro a calcular: "prima", "delta", "gamma", "theta", "vega", "rho".
 * @param {number} s Precio spot del subyacente.
 * @param {number} k Precio de ejercicio.
 * @param {number} T Tiempo hasta el vencimiento (en años).
 * @param {number} r Tasa libre de riesgo.
 * @param {number} q Tasa de dividendos o costo de carry.
 * @param {number} sigma Volatilidad del subyacente (decimal).
 * @param {string} CallPut "call" o "put".
 * @param {string} tipo_subyacente "accion" o "futuro".
 * @return Valor del parámetro solicitado.
 * @customfunction
 */
function BS93Amer(parametro, s, k, T, r, q, sigma, CallPut, tipo_subyacente) {
  parametro = parametro.toLowerCase();
  CallPut = CallPut.toLowerCase();
  tipo_subyacente = tipo_subyacente.toLowerCase();

  let b;
  if (tipo_subyacente === "accion") {
    b = r - q;
  } else if (tipo_subyacente === "futuro") {
    b = 0;
  } else {
    throw new Error("tipo_subyacente debe ser 'accion' o 'futuro'");
  }

  const h = 0.01; // paso para diferencias finitas

  function f(s, k, T, r, b, sigma) {
    return BS93AmericanOption(CallPut, s, k, T, r, b, sigma);
  }

  switch (parametro) {
    case "prima":
      return f(s, k, T, r, b, sigma);
    case "delta":
      return (f(s + h, k, T, r, b, sigma) - f(s - h, k, T, r, b, sigma)) / (2 * h);
    case "gamma":
      return (f(s + h, k, T, r, b, sigma) - 2 * f(s, k, T, r, b, sigma) + f(s - h, k, T, r, b, sigma)) / (h * h);
    case "theta":
      return (f(s, k, T - h, r, b, sigma) - f(s, k, T, r, b, sigma)) / h;
    case "vega":
      return (f(s, k, T, r, b, sigma + h) - f(s, k, T, r, b, sigma - h)) / (2 * h)/100;
    case "rho":
      return (f(s, k, T, r + h, b, sigma) - f(s, k, T, r - h, b, sigma)) / (2 * h)/100;
    default:
      return "Parámetro inválido";
  }
}

function BS93AmericanOption(optionType, S, X, T, r, b, sigma) {
  if (optionType.toLowerCase() === "put") {
    return BS93AmericanCall(X, S, T, r - b, -b, sigma);
  } else if (optionType.toLowerCase() === "call") {
    return BS93AmericanCall(S, X, T, r, b, sigma);
  } else {
    throw new Error("El tipo de opción debe ser 'call' o 'put'");
  }
}

function BS93AmericanCall(S, X, T, r, b, sigma) {
  if (b >= r) return BS_EuropeanCall(S, X, T, r, b, sigma);

  let beta = (0.5 - b / (sigma * sigma)) + Math.sqrt(Math.pow(b / (sigma * sigma) - 0.5, 2) + 2 * r / (sigma * sigma));
  let BInfinity = beta / (beta - 1) * X;
  let B0 = Math.max(X, r / (r - b) * X);
  let hT = -(b * T + 2 * sigma * Math.sqrt(T)) * (B0 / (BInfinity - B0));
  let I = B0 + (BInfinity - B0) * (1 - Math.exp(hT));
  let alpha = (I - X) * Math.pow(I, -beta);

  if (S >= I) return S - X;

  return alpha * Math.pow(S, beta)
    - alpha * phi(S, T, beta, I, I, r, b, sigma)
    + phi(S, T, 1, I, I, r, b, sigma)
    - phi(S, T, 1, X, I, r, b, sigma)
    - X * phi(S, T, 0, I, I, r, b, sigma)
    + X * phi(S, T, 0, X, I, r, b, sigma);
}

function phi(S, T, gamma, H, I, r, b, sigma) {
  let lambda = -r + gamma * b + 0.5 * gamma * (gamma - 1) * sigma * sigma;
  let d = -(Math.log(S / H) + (b + (gamma - 0.5) * sigma * sigma) * T) / (sigma * Math.sqrt(T));
  let kappa = 2 * b / (sigma * sigma) + (2 * gamma - 1);
  return Math.exp(lambda * T) * Math.pow(S, gamma) *
    (normCDF(d) - Math.pow(I / S, kappa) * normCDF(d - 2 * Math.log(I / S) / (sigma * Math.sqrt(T))));
}

function BS_EuropeanCall(S, X, T, r, b, sigma) {
  let d1 = (Math.log(S / X) + (b + sigma * sigma / 2) * T) / (sigma * Math.sqrt(T));
  let d2 = d1 - sigma * Math.sqrt(T);
  return S * Math.exp((b - r) * T) * normCDF(d1) - X * Math.exp(-r * T) * normCDF(d2);
}

function normCDF(x) {
  return (1.0 + erf(x / Math.sqrt(2))) / 2.0;
}

function erf(x) {
  let a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  let a4 = -1.453152027, a5 = 1.061405429;
  let p = 0.3275911;
  let sign = x < 0 ? -1 : 1;
  x = Math.abs(x);
  let t = 1 / (1 + p * x);
  let y = 1 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * Math.exp(-x * x);
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


