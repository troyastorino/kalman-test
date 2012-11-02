(ns kalman-sim.kalman-filter
  (:use [incanter.core]
        [incanter.stats]))

(defprotocol Kalman-Filter
  (g [u x] "Models the state transitions\nu: command vector\nx: state vector")
  (G [u x] "Generates the Jacobian matrix of g\nu: command vector\nx: state vector")
  (h [x] "Models measurment\nx: state vector")
  (H [x] "Generates the Jacobian matrix of h\n x: state vector"))

(defn generate-EKF
  "Returns an extended kalman filter function with params [state cov command meas]\n
  state-fn: a function modeling state transition. Params must be [command state]\n
  meas-fn: a function modeling measurment. Params must be [state]\n
  state-jacobian: a function returning the Jacobian matrix of the state function. Params must be [command state]\n
  meas-jacobian: a function returning the Jacobian matrix of the meas function. Params must be [state]\n
  state-cov: a covariance matrix for the state transition function\n
  meas-cov: a covariance matrix for the measurment function"
  [state-fn meas-fn state-jacobian meas-jacobian state-cov meas-cov]
  (fn [x cov command meas]
    (let [x-pred (state-fn command x) ; predicted state estimate
          G (state-jacobian command x) ; state transition jacobian
          cov-pred (plus (mmult G cov (trans G)) state-cov) ; predicted state estimate covariance
          H (meas-jacobian x-pred) ; measurment jacobian
          KG (mmult cov-pred (trans H) (solve (plus (mmult H cov-pred (trans H)) meas-cov))) ; kalman gain
          meas-pred (meas-fn x-pred) ; predicted measurment
          new-x (plus x-pred (mmult KG (minus meas meas-pred))) ; updated state estimate
          new-cov (mmult (minus (identity-matrix (count KG)) (mmult KG H)) cov-pred)] ; updated estimate covariance
      [new-x new-cov])))