def test_pedal():
    # Test with a vertical line
    assert pedal((0, 0), (0, 10), (5, 5)) == (0, 5), "Failed on vertical line"
    # Test with a horizontal line
    assert pedal((0, 0), (10, 0), (5, 5)) == (5, 0), "Failed on horizontal line"
    # Test with a diagonal line
    assert pedal((0, 0), (10, 10), (5, 0)) == (2.5, 2.5), "Failed on diagonal line"