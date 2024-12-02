document.addEventListener("DOMContentLoaded", function () {
    const img = document.querySelector(".theme-toggle-image");

    if (!img) return; // Exit if no image is found

    // Function to update the image based on the current theme
    function updateImage() {
        const theme = document.body.getAttribute("data-md-color-scheme");
        img.src = theme === "slate" ? "images/NPLinker_standard_white.svg" : "images/NPLinker_standard_black.svg";
    }

    // Initial update
    updateImage();

    // Observe changes to the `data-md-color-scheme` attribute
    const observer = new MutationObserver(updateImage);
    observer.observe(document.body, {
        attributes: true,
        attributeFilter: ["data-md-color-scheme"], // Monitor changes to the theme attribute
    });
});