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



