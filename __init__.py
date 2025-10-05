from .bloom_watch import BloomWatchPlugin


def classFactory(iface) -> BloomWatchPlugin:
    return BloomWatchPlugin(iface)
