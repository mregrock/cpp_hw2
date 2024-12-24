# benchmark

After replacing std::ranges::find with simple if-else, the performance increased by 100% on 200 ticks benchmark.

Run 1:
  Total time: 6649.000 ms
  Average time per tick: 33.245 ms
  Ticks per second: 30.080

Run 2:
  Total time: 6521.000 ms
  Average time per tick: 32.605 ms
  Ticks per second: 30.670

Run 3:
  Total time: 6530.000 ms
  Average time per tick: 32.650 ms
  Ticks per second: 30.628

Run 4:
  Total time: 6530.000 ms
  Average time per tick: 32.650 ms
  Ticks per second: 30.628

Run 5:
  Total time: 6523.000 ms
  Average time per tick: 32.615 ms
  Ticks per second: 30.661

Average across 5 runs:
  Total time: 6550.600 ms
  Average time per tick: 32.753 ms
  Ticks per second: 30.533



