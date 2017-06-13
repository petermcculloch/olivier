(ns olivier.math
  (:require [clojure.math.numeric-tower :as math]))

(defn factorization [number]
  (let [divisible? (fn [a b] (= (rem a b) 0))]
    (loop [remaining number
           primes []
           divisor  2]
      (cond
        (= remaining 1) primes
        (divisible? remaining divisor) (recur (/ remaining divisor)
                                              (conj primes divisor)
                                              divisor)
        :else (recur remaining
                     primes
                     (inc divisor))))))

(defn clamp
  ([x] (clamp x 0 1))
  ([x x-min x-max] (min (max x-min x) x-max)))

(defn wrap
  "Wraps a number.  Inclusive? defaults to true."
  ([x] (wrap x 0 1 true))
  ([x xmin xmax] (wrap x xmin xmax true))
  ([x xmin xmax inclusive?]
   (let [offset (if inclusive? 1 0)
         [xmin xmax] [(min xmin xmax) (+ (max xmin xmax) offset)]]
     (+ (mod (- x xmin) (- xmax xmin)) xmin))))

(defn fold
  "Folds a number.  Inclusive? defaults to true."
  ([x] (wrap x 0 1 true))
  ([x xmin xmax] (wrap x xmin xmax true))
  ([x xmin xmax inclusive?]
   (let [offset (if inclusive? 0 1)
         [xmin xmax] [(min xmin xmax) (max xmin xmax)]
         dist (max 1 (- xmax (+ xmin offset)))
         dbl-dist (* 2 dist)
         x (mod (- x xmin) dbl-dist)]
     (+ (if (> x dist) (- dbl-dist x) x) xmin))))

(defn fold-dbl
  "Folds a number, doubling the endpoints.  Inclusive? defaults to true."
  ([x] (wrap x 0 1 true))
  ([x xmin xmax] (wrap x xmin xmax true))
  ([x xmin xmax inclusive?]
   (let [offset (if inclusive? 1 0)
         [xmin xmax] [(min xmin xmax) (max xmin xmax)]
         dist (max 1 (+ (- xmax xmin) offset))
         dbl-dist (* 2 dist)
         x (mod (- x xmin) dbl-dist)]
     (println x xmin xmax dist dbl-dist)
     (+ (if (> x (dec dist)) (- (dec dbl-dist) x) x) xmin))))

(defn numer [x] (if (ratio? x) (numerator x) x))
(defn denom [x] (if (ratio? x) (denominator x) 1))

(defn find-fraction
  "Returns fractional approximation with the lowest denominator that meets the specified epsilon.
   If epsilon is not met, incEpsilonFn is applied to epsilon until epsilon is met.
   If not incEpsilonFn is supplied, epsilon will be doubled."
  ([x] (find-fraction x 1 100 0.00001 #(* 2 %)))
  ([x upper-limit] (find-fraction x 1 upper-limit 0.00001 #(* 2. %)))
  ([x lower-limit upper-limit] (find-fraction x lower-limit upper-limit 0.00001 #(* 2 %)))
  ([x lower-limit upper-limit epsilon] (find-fraction x lower-limit upper-limit epsilon #(* 2 %)))
  ([x lower-limit upper-limit epsilon inc-epsilon-fn]
   (if (not (float? x))
     x
     ; Calculate this once
     (let [v (mapv #(vector % (Math/abs (- (Math/round (* x %)) (* x %)))) (range lower-limit (inc upper-limit)))]
       (loop [eps epsilon]
         (let [c (ffirst (filterv #(< (second %) eps) v))]
           (if (not (nil? c))
             (/ (int (Math/round (* x c))) c)
             (recur (inc-epsilon-fn eps)))))))))

(defn- make-in-range [p1 p2]
  #(and (>= % (:x p1)) (<= % (:x p2))))

(defn- make-linear-interpolation
  "Linear interpolation between x1 y1 and x2 y2"
  [p1 p2]
  (fn [x] (let [x1 (:x p1)
                y1 (:y p1)
                x2 (:x p2)
                y2 (:y p2)
                dist (/ (- x2 x) (- x2 x1))]
            (+ (* y2 (- 1 dist)) (* y1 dist)))))

(defn make-linear-bpf [m bounds-mode]
  "M is a vector of points"
  (let [pairs (partition 2 1 m)
        fns (map #(apply make-linear-interpolation %) pairs)
        lower (first m)
        upper (last m)
        in-ranges (map #(apply make-in-range %) pairs)
        make-x (fn [x] ((cond (= :clip bounds-mode) clamp
                              (= :wrap bounds-mode) wrap)
                         x (:x lower) (:x upper)))
        fn-ranges (map vector in-ranges fns)
        ]
    ;[fn-ranges ]))

    (fn [x]
      (let [x (make-x x)]
        ((second (first (drop-while #(not ((first %) x)) fn-ranges))) x)))))


