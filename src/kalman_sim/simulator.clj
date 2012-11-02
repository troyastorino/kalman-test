(ns kalman-sim.simulator
  (:use [incanter.core]
        [incanter.stats]
        [kalman-sim.kalman-filter]))

(def delta-t 1)
(def max-measurments 1e2)

(def state (ref []))

(defrecord Log-Entry [mu sigma state])

(def history (ref []))

(defn add-log-entry! [mean sigma state]
  (dosync
   (alter history conj (Log-Entry. mean sigma state))))

(defn reset-log! []
  (dosync
   (ref-set history [])))

(defn set-state! [x]
  (dosync
   (ref-set state x)))

(defn gen-update-state
  "Returns a funtion that generates the next state that takes arguments [u x]\n
  g: state transition equation\nQcovariance in state transition"
  [g R]
  (fn [u x]
    (let [mu (g u x)
          v (sample-mvn 1 :mean mu :sigma R)]
      v)))

(defn gen-sample
  "Returns a funtion that generates a sample measurment which takes arguments [x]\n
  h: measurment funciton\nQ: covariance in measurments"
  [h Q]
  (fn [x]
    (let [mu (h x)
          v (sample-mvn 1 :mean mu :sigma Q)]
      v)))

(defn run-simulation! [g h G H R Q mu0 sig0 planner init-state]
  (reset-log!)
  (set-state! init-state)
  (let [EKF (generate-EKF g h G H R Q)
        sample (gen-sample h Q)
        update-state (gen-update-state g R)]
    (add-log-entry! mu0 sig0 @state)
    (loop [mu mu0 sig sig0 i 0]
      (let [x @state
            command (planner mu sig)
            x-new (update-state command x)]
        (set-state! x-new)
        (let [meas (sample x-new)
              [mu-new sig-new] (EKF mu sig command meas)]
        (add-log-entry! mu-new sig-new x-new)
        (if (and (< i max-measurments) (> (trace sig-new) 0.05))
          (do
            (recur mu-new sig-new (inc i)))))))))