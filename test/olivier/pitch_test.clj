(ns olivier.pitch-test
  (:require [clojure.test :as test])
  (:use [olivier.pitch]))

(def algochords
  ; Chords from David Cope's "The Algorithmic Composer"
  {
   :3-20   {:c1 [60 63 65 68]                               ; Chords 1 & 2 from fig 3.20,  p 89
            :c2 [62 65 63 67]                               ;
            :r  [[[0 0] 2 0]                                ; Rule
                 [[3 3] 2 1]
                 [[5 1] -2 2]
                 [[8 5] -1 3]]}


   :3-21-a {:c1 [36 50 69 77],                               ; Chords 1 & 2 from fig 3.21a, p 90
            :c2 [38 48 71 76],                               ;
            :r  [[[3 3] 2  0]
                 [[5 1] -2 1]
                 [[0 0] 2  2]
                 [[8 5] -1 3]]}

   :3-21-c {:c1 [45 61 64 78],                               ;  ---              fig 3.21c, p 90
            :c2 [50 58 67 76],                               ;
            :r  [[[8 0] 5  0]
                 [[0 8] -3 1]
                 [[3 5] 3  2]
                 [[5 2] -2 3]]}

   :3-21-e {:c1 [36 50 69 77],                              ;  ---              fig 3.21e, p 90
            :c2 [39 48 66 82],                              ;
            :r  [[[3 5] 3 0]
                 [[5 2] -2 1]
                 [[0 8] -3 2]
                 [[8 0] 5 3]]}
   })

(test/deftest create-rule-test
  (letfn [(make-rule-test
            [m]
            (let [{:keys [c1 c2 r]} m]
              (= (create-rule c1 c2) r)))]
    (test/testing "Create rules from chord transition"
      (map #(test/is (make-rule-test %)) (keys algochords)))))


(test/deftest test-apply-chord-transition
  (test/testing "Test chord transition")
  (letfn [(chord-trans [m]
            (let [{:keys [c1 c2 r]} m]
              (= (apply-chord-transition c1 r) c2)))]
    (test/testing "Create rules from chord transition"
      (map #(test/is (chord-trans %)) (vals algochords)))))

