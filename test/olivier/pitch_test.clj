(ns olivier.pitch-test
  (:require [clojure.test :as test])
  (:use [olivier.pitch]))

(def algochords
  ; Chords from David Cope's "The Algorithmic Composer"
  {:3-20-1      [60 63 65 68] ; Chords 1 & 2 from fig 3.20,  p 89
   :3-20-2      [62 65 63 67] ;
   :3-20-rule   [[[0 0] 2 0]
                 [[3 3] 2 1]
                 [[5 1] -2 2]
                 [[8 5] -1 3]]

   :3-21-a1     [36 50 69 77] ; Chords 1 & 2 from fig 3.21a, p 90
   :3-21-a2     [38 48 71 76] ;
   :3-21-a-rule [[[3 3] 2 0]
                 [[5 1] -2 1]
                 [[0 0] 2 2]
                 [[8 5] -1 3]]

   :3-21-c1     [45 61 64 78] ;  ---              fig 3.21c, p 90
   :3-21-c2     [50 58 67 76] ;
   :3-21-c-rule [[[8 0] 5 0]
                 [[0 8] -3 1]
                 [[3 5] 3 2]
                 [[5 2] -2 3]]

   :3-21-e1     [36 50 69 77] ;  ---              fig 3.21e, p 90
   :3-21-e2     [39 48 66 82] ;
   :3-21-e-rule [[[3 5] 3 0]
                 [[5 2] -2 1]
                 [[0 8] -3 2]
                 [[8 0] 5 3]]})

(defn- make-rule-test [k1 k2 k3]
  (= (create-rule (algochords k1) (algochords k2))
     (algochords k3)))

(test/deftest create-rule-test
         (test/testing "Create rules from chord transition"
                  (test/is (make-rule-test :3-20-1  :3-20-2  :3-20-rule))
                  (test/is (make-rule-test :3-21-a1 :3-21-a2 :3-21-a-rule))
                  (test/is (make-rule-test :3-21-c1 :3-21-c2 :3-21-c-rule))
                  (test/is (make-rule-test :3-21-e1 :3-21-e2 :3-21-e-rule))))
