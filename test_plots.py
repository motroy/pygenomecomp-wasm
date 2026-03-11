import asyncio
from playwright.async_api import async_playwright
import time

async def main():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        print("Navigating to local server...")
        await page.goto("http://localhost:8000/")

        # Click the 'Load Real-World Example' button
        print("Clicking Load Real-World Example...")
        await page.click("#btn-load-ndm")

        # We need to click the visual switch label because the input itself is screen-reader-only
        print("Clicking 'Show gene names' label...")
        await page.evaluate('document.getElementById("show-gene-names").click()')

        # Check 'Resistance genes only'
        print("Clicking 'Resistance genes only' label...")
        await page.evaluate('document.getElementById("show-res-genes").click()')

        # Run Comparison
        print("Running comparison...")
        await page.click("#btn-run")

        # Wait for SVG container to populate
        print("Waiting for results to render...")
        await page.wait_for_selector("#svg-container svg", timeout=120000)

        print("Taking screenshot of circular view...")
        await page.screenshot(path="plot_circular.png", full_page=True)

        # Click linear view
        print("Switching to linear view...")
        await page.click("#btn-view-linear")
        # wait a bit
        await asyncio.sleep(2)
        print("Taking screenshot of linear view...")
        await page.screenshot(path="plot_linear.png", full_page=True)

        # Click alignment view
        print("Switching to alignment view...")
        await page.click("#btn-view-alignment")
        # wait a bit
        await asyncio.sleep(2)
        print("Taking screenshot of alignment view...")
        await page.screenshot(path="plot_alignment.png", full_page=True)

        print("Done!")
        await browser.close()

asyncio.run(main())
