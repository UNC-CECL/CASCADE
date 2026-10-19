# first line: 327
def _fetch_tile(tile_url, wait, max_retries, headers: dict[str, str]):
    array = _retryer(tile_url, wait, max_retries, headers)
    return array
