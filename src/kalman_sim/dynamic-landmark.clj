(ns kalman-sim.dynamic-landmark
  (:use [incanter.core]
        [incanter.stats]
        [incanter.charts]
        [kalman-sim.simulator]))

(def size 100)

(defn init-state []
  (let [xr (- (rand (* 2 size)) size)
        yr (- (rand (* 2 size)) size)
        theta (- (rand (* 2 Math/PI)) Math/PI)
        xl (- (rand (* 2 size)) size)
        yl (- (rand (* 2 size)) size)]
    (dosync
     (ref-set state [xr yr theta xl yl]))))

(def Q (matrix [[1e-3 0 0 0 0] 
                [0 1e-3 0 0 0]
                [0 0 1e-4 0 0]
                [0 0 0 1e-3 0]
                [0 0 0 0 1e-3]])) ; Cov in executing commands
(def R (matrix [[100 0]
                [0 50]])) ; Cov in measurment

(defn wrap-angle [theta]
  (if (> theta Math/PI)
    (- theta (* 2 Math/PI))
    (if (<= theta (- Math/PI) )
      (+ theta (* 2 Math/PI))
      theta)))

(defn model-g
  "The function modeling state transitions\n
  u: command vector\n
  x: state vector"
  [u x]
  (let [[xr yr theta xl yl] x
        [v omega] u
        nxr (+ xr  (* v (Math/cos theta) delta-t))
        nyr (+ yr  (* v (Math/sin theta) delta-t))
        ntheta (wrap-angle (+ theta (* omega delta-t)))]
    [nxr nyr ntheta xl yl]))

(defn true-g [u x]
  (let [[xr yr theta xl yl] (model-g u x)]
    [xr yr theta (dec xl) (dec yl)]))

(defn G
  "The Jacobian of g\n
  u: command vector\n
  x: state vector"
  [u x]
  (let [[x y theta] x
        [v omega] u]
    (matrix [[1 0 (* -1 v (Math/sin theta) delta-t) 0 0]
             [0 1 (* v (Math/cos theta) delta-t) 0 0]
             [0 0 1 0 0]
             [0 0 0 1 0]
             [0 0 0 0 1]])))

(defn h
  "The function modeling measurment\n
  x: state vector"
  [x]
  (let [[xr yr theta xl yl] x
        xd (- xr xl)
        yd (- yr yl)
        d (Math/sqrt (+ (* xd xd) (* yd yd)))
        phi (+ Math/PI (Math/atan2 yd xd) (- theta))]
    [d phi]))

(defn H
  "The Jacobian of h\n
  x: state vector"
  [state]  (let [[xr yr theta xl yl] state
                 xd (- xr xl)
                 yd (- yr yl)
                 l2 (+ (* xd xd) (* yd yd))
                 d (Math/sqrt l2)]
             (matrix [[(/ xd d) (/ yd d) 0 (- (/ xd d)) (- (/ yd d))]
                      [(- (/ yd l2)) (/ xd l2) -1  (/ yd l2)  (- (/ xd l2))]])))

(defn planner [mu sigma]
  [0 0])

(defn dist [[x y theta xl yl]]
  (Math/sqrt (+ (* x x) (* y y))))

(defn sample-data [log-entry]
  (sample-mvn 100 :mean (:mu log-entry) :sigma (:sigma log-entry)))

(defn plot-belief []
  (let [entries @history
        plot (box-plot (map dist (sample-data (first entries)))
                       :x-label "Number of Measurments"
                       :y-label "Distance from Landmark"
                       :title "Belief")] 
    (doseq [entry (rest entries)]
      (add-box-plot plot (map dist (sample-data entry))))
    (view plot)))

(defn plot-true-state []
  (let [entries @history
        plot (xy-plot (range (count entries)) (map #(dist (:state %1)) entries)
                   :x-label "Number of Measurments"
                   :y-label "Distance from Landmark"
                   :title "True State")]
    (set-x-label plot "Number of Measurments")
    (set-y-label plot "Distance from Landmark")
    (set-title plot "True State" )
    (view plot)))

(defn -main
  "Runs the simulation, outputs a graph when the trace of the
  covariance matrix is small enough."
  [& arg]
  (run-simulation! model-g h G H Q R [(/ size 2) (/ size 2) 0 0 0] (matrix [[size 0 0 0 0]
                                                                            [0 size 0 0 0]
                                                                            [0 0 Math/PI 0 0]
                                                                            [0 0 0 size 0]
                                                                            [0 0 0 0 size]])
                   planner (init-state))
  (plot-belief)
  (plot-true-state))